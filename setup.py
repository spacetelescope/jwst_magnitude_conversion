from setuptools import setup, find_packages
import glob

setup(name="jwst_magnitude_conversion",
      version="0.1.0",
      description="Code to estimate JWST imaging filter magnitudes on the basis of external colour and magnitude information.",
    author="Kevin Volk", author_email="volk@stsci.edu",
    keywords=["jwst", "commissioning", "niriss"],
    classifiers=['Programming Language :: Python', 'Programming Language :: Python :: 3',
                 'Development Status :: 1 - Planning', 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering :: Astronomy',
                 'Topic :: Scientific/Engineering :: Physics',
                 'Topic :: Software Development :: Libraries :: Python Modules'],
    py_modules=[x.split(".py")[0] for x in glob.glob("*py") if "setup.py" not in x],
    packages=find_packages(),

    # dependencies should be taken care of by the environment file
    install_requires=["setuptools", "matplotlib", "astropy", "configobj", "numpy"
                    ]
      )

