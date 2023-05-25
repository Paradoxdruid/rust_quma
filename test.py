import time

import python_quma
import rust_quma

start = time.perf_counter()
py_quma_output = python_quma.Quma(
    ">query\nATCGTAGTCGA", ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"
)
end = time.perf_counter()

print(py_quma_output.data)
print(f"Python version took {end - start:0.6f} seconds")

start = time.perf_counter()
rust_quma_output = rust_quma.Quma(
    ">query\nATCGTAGTCGA", ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"
)
end = time.perf_counter()

print(rust_quma_output.data)
print(f"Rust version took {end - start:0.6f} seconds")
