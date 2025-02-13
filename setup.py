from setuptools import find_packages, setup

setup(
    name="PauliPropagationQiskit",
    version="0.1.0",
    description="Small wrapper to expose PauliPropagation to Qiskit",
    author="Marc",
    author_email="marcsanzdrudis@outlook.com",
    packages=find_packages(),
    install_requires=[
        "juliapkg",
        "juliacall>=0.9.0",
        "qiskit",
        "qiskit_algorithms",
    ],
)
