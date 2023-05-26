import time

from maturin import import_hook

import python_quma

import_hook.install(release=True)

import rust_quma  # noqa # type: ignore

start = time.perf_counter()
for _ in range(100):
    __ = python_quma.Quma(
        ">query\nATCGTAGTCGA", ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"
    )

py_quma_output = python_quma.Quma(
    ">query\nATCGTAGTCGA", ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"
)

end = time.perf_counter()

pytime = end - start

print(py_quma_output.values)
print(f"Python version took {pytime:0.6f} seconds")

start = time.perf_counter()
for _ in range(100):
    _ = rust_quma_output = rust_quma.Quma(
        ">query\nATCGTAGTCGA", ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"
    )

rust_quma_output = rust_quma.Quma(
    ">query\nATCGTAGTCGA", ">query1\nATCGTAGTCGA\n>query2\nATCGATAGCATT"
)
end = time.perf_counter()

rust_time = end - start

print(rust_quma_output.values)
print(f"Rust version took {rust_time:0.6f} seconds")

print(f"Rust version is {pytime/rust_time:0.2f} times faster")
