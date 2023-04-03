import setuptools

setuptools.setup(
    name="SDMBench", # Replace with your own username
    version="1.0.0",
    author="fangyuan zhao",
    author_email="zhaofangyuan20s@ict.ac.cn",
    description="SDMBench metrics function for python",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    install_requires = [
    'scanpy',
    'pandas',
    'numpy',
    'squidpy',
    'scikit-learn',
    'requests',
    'urllib3>=1.26.12',
    'tqdm',
    'anndata'
    ] 
)
