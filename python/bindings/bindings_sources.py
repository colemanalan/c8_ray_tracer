import os


def get_sources():
    source_dir = os.path.dirname(__file__)

    sources = [
        "cash_karp_bindings.cpp",
        "environment_bindings.cpp",
        "logger_bindings.cpp",
        "path_bindings.cpp",
        "plane_bindings.cpp",
        "pywrapper.cpp",
        "ray_tracing_table_bindings.cpp",
        "raytracer_bindings.cpp",
        "vec_bindings.cpp",
    ]

    return [os.path.join(source_dir, src) for src in sources]


if __name__ == "__main__":
    print(" ".join(get_sources()))
