import os


def get_sources():
    source_dir = os.path.dirname(__file__)

    sources = [
        "pywrapper.cpp",
        "logger_bindings.cpp",
        "tracer_bindings.cpp",
        "vec_bindings.cpp",
    ]

    return [os.path.join(source_dir, src) for src in sources]


if __name__ == "__main__":
    print(" ".join(get_sources()))
