"""Memory-guarded dense Cholesky experiments for macOS; requires NumPy, psutil, clang++.

Allocate each full dense matrix, including padding. Skip zero contributions but
preserve the order of all effective operations, with separate multiply/subtract.
"""
import argparse
import ast
import ctypes
import gc
import json
import math
import os
from pathlib import Path
import re
import resource
import subprocess
import tempfile
import time
from fractions import Fraction

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
import numpy as np
import psutil

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
GiB=2**30

def load_reference():
    ns={}
    exec("import numpy as np\nimport math\nfrom fractions import Fraction",ns)
    fence=chr(96)*3
    pattern=fence+r"(?:\{code-cell\} ipython3|python)\n(.*?)"+fence
    for code in re.findall(pattern,(ROOT/"content/cholesky.md").read_text(),re.S):
        for node in ast.parse(code).body:
            if isinstance(node,ast.FunctionDef) and node.name in {"cholesky_in_place","cholesky_breakdown_matrix"}:
                exec(compile(ast.Module(body=[node],type_ignores=[]),"content/cholesky.md","exec"),ns)
    return ns

def verify_input_value_classes(q, dtype):
    """Every input entry belongs to one of these exact rational value classes."""
    p=np.finfo(dtype).nmant+1
    k=4**q
    eta=Fraction(12)*Fraction(2)**(q-p)
    delta=k*eta
    t=Fraction(3)*Fraction(2)**(-(p+q+2)//2)
    root=2**q
    values=[Fraction(0),Fraction(4),2*t,Fraction(3),
            Fraction(3,root),-Fraction(3,root),eta,
            Fraction(9,4)+eta+delta,
            Fraction(9,4*root)+eta,-Fraction(9,4*root)+eta]
    assert all(Fraction(float(dtype(float(v))))==v for v in values)
    assert t*t==Fraction(9,8)*(Fraction(2)**(1-p)/root)

def swap_counters():
    # psutil's swap sysctl is restricted by the desktop sandbox; vm_stat works.
    output=subprocess.check_output(["vm_stat"],text=True)
    page_size=int(re.search(r"page size of (\d+) bytes",output).group(1))
    return {name:int(re.search(name+r":\s+(\d+)",output).group(1))*page_size
            for name in ("Swapins","Swapouts")}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=HERE/"results.json")
    parser.add_argument("--max-gib",type=float,default=3.0)
    parser.add_argument("--reserve-gib",type=float,default=2.0)
    parser.add_argument("--case",type=int,nargs=3,metavar=("BITS","Q","N"),help="Run one case in an isolated process")
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="nla-cholesky-") as tmp:
        libpath=Path(tmp)/"kernel.dylib"
        cmd=["clang++","-O3","-std=c++17","-ffp-contract=off","-fno-fast-math",
             "-dynamiclib",str(HERE/"kernel.cpp"),"-o",str(libpath)]
        subprocess.run(cmd,check=True)
        lib=ctypes.CDLL(str(libpath))
        for bits in (32,64):
            init=getattr(lib,f"init{bits}")
            chol=getattr(lib,f"chol{bits}")
            init.argtypes=[ctypes.c_void_p,ctypes.c_int,ctypes.c_int]
            init.restype=None
            chol.argtypes=[ctypes.c_void_p,ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_int64)]
            chol.restype=ctypes.c_int
        reference=load_reference()
        report={"date":time.strftime("%Y-%m-%d"),
                "total_ram_gib":psutil.virtual_memory().total/GiB,
                "numpy":np.__version__,
                "compile_flags":cmd[1:6],
                "rounding":"separate multiplication/subtraction, no FMA; zero contributions skipped",
                "memory_policy":{"max_matrix_gib":args.max_gib,"reserve_available_gib":args.reserve_gib},
                "runs":[],"skipped":[]}
        def save():
            args.output.write_text(json.dumps(report,indent=2)+"\n")
        cases=[(64,3,576),(32,4,2304),(64,5,9216),
               (64,5,12288),(64,5,14336),(64,5,16384),
               (64,5,17408),(64,5,18432),(64,5,19456),
               (64,5,20480),(64,5,23040),(64,7,147456)]
        if args.case:
            cases=[tuple(args.case)]
        for bits,q,n in cases:
            dtype=np.dtype(f"float{bits}")
            p=53 if bits==64 else 24
            assert q>=3 and (p+q)%2==0 and 3*q+8<=p
            k=4**q
            assert n >= 9*k
            verify_input_value_classes(q,dtype.type)
            size=n*n*dtype.itemsize
            available=psutil.virtual_memory().available
            extra=size if n<=2304 else 0
            needed=size+extra+128*2**20
            if size>args.max_gib*GiB or needed>available-args.reserve_gib*GiB:
                item={"bits":bits,"q":q,"n":n,"matrix_gib":size/GiB,
                      "available_gib":available/GiB,"reason":"memory guard"}
                report["skipped"].append(item)
                print("SKIP",json.dumps(item),flush=True)
                save()
                continue
            print(f"START float{bits} q={q} n={n}: matrix={size/GiB:.3f} GiB, available={available/GiB:.3f} GiB",flush=True)
            before_swap=swap_counters()
            start=time.perf_counter()
            a=np.empty((n,n),dtype=dtype,order="F")
            getattr(lib,f"init{bits}")(a.ctypes.data,n,q)
            build_s=time.perf_counter()-start
            ref=None
            if n<=2304:
                ref=reference["cholesky_breakdown_matrix"](q,dtype.type)
                assert a.tobytes(order="F")==ref.tobytes(order="F"),"constructor differs from notes"
            bad=ctypes.c_int64()
            start=time.perf_counter()
            step=getattr(lib,f"chol{bits}")(a.ctypes.data,n,q,ctypes.byref(bad))
            factor_s=time.perf_counter()-start
            assert step>0 and bad.value==0,(step,bad.value)
            pivot=float(a[step-1,step-1])
            if ref is not None:
                try:
                    reference["cholesky_in_place"](ref)
                except np.linalg.LinAlgError as error:
                    assert str(step) in str(error),(step,error)
                else:
                    raise AssertionError("reference unexpectedly succeeded")
                assert a.tobytes(order="F")==ref.tobytes(order="F"),"native result differs bitwise from notes"
                del ref
            xf=Fraction(2)**(3*q-p)
            bound=(Fraction(17,2)+39*xf)*(1+Fraction(17,72)/xf)
            after_swap=swap_counters()
            result={"bits":bits,"q":q,"base_n":9*k,"n":n,"padded":n!=9*k,
                    "matrix_gib":size/GiB,"x":float(xf),
                    "condition_upper_bound_exact":str(bound),
                    "condition_upper_bound_float":float(bound),
                    "failure_step":step,"pivot":pivot,"pivot_hex":pivot.hex(),
                    "schur_entry_mismatches":bad.value,
                    "input_value_classes_exact":True,
                    "reference_bitwise_match":True if n<=2304 else None,
                    "build_seconds":build_s,"factor_seconds":factor_s,
                    "process_peak_rss_gib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/GiB,
                    "available_after_gib":psutil.virtual_memory().available/GiB,
                    "system_swapin_delta_bytes":after_swap["Swapins"]-before_swap["Swapins"],
                    "system_swapout_delta_bytes":after_swap["Swapouts"]-before_swap["Swapouts"]}
            report["runs"].append(result)
            print("RESULT",json.dumps(result),flush=True)
            save()
            del a
            gc.collect()
        save()
if __name__=="__main__":
    main()
