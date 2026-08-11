<div align="center">
    <a name="readme-top"></a>
    <h1>
        PyMMseqs 🚀
    </h1>


![GitHub Actions](https://img.shields.io/github/actions/workflow/status/peymanvahidi/pymmseqs/release.yaml?style=plastic&logo=github-actions&label=CI)
[![License](https://img.shields.io/github/license/peymanvahidi/pymmseqs?style=plastic&color=orange&logo=github&label=License)](./LICENSE)
[![PyPI Version](https://img.shields.io/pypi/v/pymmseqs?style=plastic&color=4b8bbe&logo=pypi&logoColor=white&label=PyPI)](https://pypi.org/project/pymmseqs/)
[![Docker](https://img.shields.io/badge/Docker-GHCR-0db7ed?style=plastic&logo=docker&label=Docker)](https://github.com/peymanvahidi/pymmseqs/pkgs/container/pymmseqs)

PyMMseqs is a powerful Python wrapper for [MMseqs2](https://github.com/soedinglab/MMseqs2). It seamlessly integrates MMseqs2’s advanced functionality into your Python workflows, allowing you to effortlessly execute MMseqs2 commands and parse their outputs into convenient Python objects for further analysis. Whether you're clustering sequences, searching databases, or analyzing large-scale biological data, PyMMseqs simplifies the process while maintaining the performance and flexibility of MMseqs2.

</div>

---

## 🎯 5-Minute Tour of PyMMseqs

Want to see the power of PyMMseqs in action?

Try our quick tour on Google Colab and experience the full power of PyMMseqs, all in just 5 minutes! 

**Ready to dive in?**

[![PyMMseqs Quick Tour](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1xR78P1OZqBDKTM5feTzJG_MHn-TlOdNj?usp=sharing)

---

## 🗝️ Features

- **Seamless Integration**: Execute MMseqs2 commands directly within your Python code, eliminating the need for shell scripting or external command-line tools.
- **Output Parsing**: Convert MMseqs2 outputs into Python objects (e.g., Pandas DataFrame, generators, dictionaries) for easy manipulation and analysis.
- **High Performance**: Leverage the speed and efficiency of MMseqs2 while enjoying the flexibility of Python.
- **Cross-Platform**: Use PyMMseqs via pip or Docker, ensuring compatibility across different environments. PyMMseqs works seamlessly on Linux and macOS.

> **Note**: Windows users should either use Windows Subsystem for Linux (WSL) or Docker to run PyMMseqs.

---

## 🛠️ Installation

PyMMseqs can be installed in two ways: via pip (recommended for most users) or using a Docker image (ideal for reproducible environments).

### Installing via pip
The `pymmseqs` package is currently available on PyPI. To install it, use the following command:

```bash
pip install pymmseqs
```

### Using Docker Image
For users who prefer not to install PyMMseqs locally or want a pre-configured environment, a Docker image is available on GitHub Container Registry (GHCR).

#### Debian-based Image
To pull the Debian-based Docker image, run:

```bash
docker pull ghcr.io/heispv/pymmseqs:latest-debian
```
> **Note**: If you want to use a specific version of PyMMseqs, you can replace `latest` with the desired version.

> **Tip**: Using Docker ensures that all dependencies, including MMseqs2, are pre-installed and configured, making it ideal for reproducible workflows.

---

## 🚀 Example Usage

Here's a simple example to get you started with PyMMseqs. This example demonstrates how to perform sequence clustering and parse the results.

If you were using MMseqs2 directly in the terminal, you would run the following command to cluster sequences:

```bash
mmseqs easy-cluster human.fasta human_clust tmp --min-seq-id 0.9
```

With PyMMseqs, you can achieve the same result directly in Python, and parse the output to Python objects for further analysis.

```python
from pymmseqs.commands import easy_cluster

# Perform clustering on a FASTA file (equivalent to the terminal command above)
human_cluster = easy_cluster("human.fasta", "human_clust", "tmp", min_seq_id=0.9)

# Get results as a Python generator for easy processing
cluster_gen = human_cluster.to_gen()

# Let's get the representative sequence of a cluster with more than 100 members
for cluster in cluster_gen:
    if len(cluster["members"]) > 100:
        print(f"Representative sequence of a large cluster: {cluster['rep']}")
        break
```

---

## 📖 Documentation

For detailed usage instructions, advanced examples, and API references, please visit the [PyMMseqs Wiki](https://github.com/heispv/pymmseqs/wiki).

---

## 🔧 Prerequisites

To use PyMMseqs, you only need:
- **Python**: Version 3.10 or higher.

> **Note**: All other dependencies, including MMseqs2, are automatically installed when you install `pymmseqs` via pip or use the Docker image.

### Where the MMseqs2 binary comes from

Wheels for Linux (x86_64, aarch64) and macOS (Apple Silicon, Intel) ship the MMseqs2 binary
inside the package, so `pip install pymmseqs` needs no extra step. PyMMseqs looks for it in
this order:

1. `MMSEQS2_PATH`, if set — use this to point at your own MMseqs2 build.
2. The binary bundled in the installed package.
3. A copy previously downloaded to `~/.cache/pymmseqs/` (respects `XDG_CACHE_HOME`).
4. Otherwise it downloads the pinned MMseqs2 release into that cache on first use.

Step 4 only happens on platforms we ship no wheel for. If your machine has no internet access,
install MMseqs2 separately and set `MMSEQS2_PATH`.

> **⚠️ Linux x86_64 requires AVX2** (Intel Haswell / AMD Excavator, 2013 or newer). The bundled
> x86_64 binary is the AVX2 build, and a Python wheel cannot express a CPU-feature requirement.
> On an older CPU it fails with `Illegal instruction` (SIGILL). If that happens, grab the
> `sse41` build from the [MMseqs2 releases](https://github.com/soedinglab/MMseqs2/releases)
> and set `MMSEQS2_PATH` to it.

---

## 🤝 Contributing

We'd love your contributions to PyMMseqs! Simply fork, branch, commit, push, and open a PR.

For bug reports, feature requests, or questions, please open an issue on the [GitHub Issues page](https://github.com/heispv/pymmseqs/issues).

---

## 📜 License

PyMMseqs is licensed under the [MIT License](LICENSE).

---

## 🌟 Support

If you find PyMMseqs useful, please consider giving the repository a star on GitHub! ⭐

It helps others discover the project and motivates further development.

For questions, feedback, or support, feel free to open an issue or contact the maintainers.

<p align="right" style="font-size: 14px; color: #555; margin-top: 20px;">
    <a href="#readme-top" style="text-decoration: none; color: #007bff; font-weight: bold;">
        ↑ Back to Top
    </a>
</p>
