import time

import rust_quma

import python_quma

start = time.perf_counter()
rust_quma_output = rust_quma.Quma(
    ">query\nATCGTAGTCGA", ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"
)
end = time.perf_counter()

print(f"Rust version took {end - start:0.4f} seconds")

start = time.perf_counter()
py_quma_output = python_quma.Quma(
    ">query\nATCGTAGTCGA", ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"
)
end = time.perf_counter()

print(f"Python version took {end - start:0.4f} seconds")