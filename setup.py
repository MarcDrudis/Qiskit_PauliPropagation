from setuptools import find_packages, setup
from setuptools.command.develop import develop
from setuptools.command.install import install


def install_pauli_propagation():
    try:
        import juliapkg

        juliapkg.add(
            "PauliPropagation",
            "293282d5-3c99-4fb6-92d0-fd3280a19750",
            url="https://github.com/MSRudolph/PauliPropagation.jl.git",
            rev="dev",
        )
        print("Julia package installed successfully.")
    except ImportError:
        print(
            "Error: juliacall could not be imported. Ensure it is installed correctly."
        )
    except Exception as e:
        print(f"Error during Julia package installation: {e}")


class InstallPauliPropagation(install):
    """Customized install command to install a Julia package."""

    def run(self):
        super().run()
        install_pauli_propagation()


class DevelopPauliPropagation(develop):
    """Customized install command to install a Julia package."""

    def run(self):
        super().run()
        install_pauli_propagation()


setup(
    name="PauliPropagationQiskit",
    version="0.1.0",
    description="Small wrapper to expose PauliPropagation to Qiskit",
    author="Marc",
    author_email="marcsanzdrudis@outlook.com",
    packages=find_packages(),
    install_requires=[
        "juliacall>=0.9.0",
        "qiskit",
    ],
    cmdclass={
        "install": InstallPauliPropagation,
        "develop": DevelopPauliPropagation,
    },
)
