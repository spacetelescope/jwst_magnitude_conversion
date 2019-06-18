from setuptools import setup, find_packages
import glob

setup(name="jwstCalibField",
      version="0.1.0",
      description="Interface to the JWST calibration field catalog.",
    author="Johannes Sahlmann", author_email="jsahlmann@stsci.edu",
    keywords=["ote", "jwst", "commissioning", "mimf", "fgs", "alignment", "niriss"],
    classifiers=['Programming Language :: Python', 'Programming Language :: Python :: 3',
                 'Development Status :: 1 - Planning', 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering :: Astronomy',
                 'Topic :: Scientific/Engineering :: Physics',
                 'Topic :: Software Development :: Libraries :: Python Modules'],
    py_modules=[x.split(".py")[0] for x in glob.glob("*py") if "setup.py" not in x],
    packages=find_packages(),

    # dependencies should be taken care of by the environment file
    install_requires=["setuptools",
                    ]
      )
