import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="Buche_Silberstein_model_2020",
    version="1.0.1",
    description="This is the Python package corresponding to: Buche, Michael R., and Meredith N. Silberstein. Statistical mechanical constitutive theory of polymer networks: The inextricable links between distribution, behavior, and ensemble. Physical Review E, 102, 012501 (2020).",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/mbuche/Buche_Silberstein_model_2020",
    author="Michael R. Buche",
    author_email="mrb289@cornell.edu",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    packages=["Buche_Silberstein_model_2020"],
    include_package_data=True,
    install_requires=["numpy", "scipy"],
)