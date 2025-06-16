# Fuzzy PSI

This project implements the Fuzzy PSI protocols presented in [Distance-Aware OT with Application to Fuzzy PSI](https://eprint.iacr.org/2025/996).

Our source code can be found inside ```sparseComp``` folder, while our benchmarks and tests can be found inside the ```tests``` folder. Additionaly, we also provide the build scripts ```build-bench.sh```,```build-tests.sh``` which can used to build our benchmarks and tests, respectively. Both of these scripts don't take any parameters. We also provide a script ```install-dependencies.sh``` that installs all the dependencies required by our project. Finally, we also provide a Dockerfile that contains all the required project dependencies and our source code. We suggest anybody simply desiring to run our benchmarks to run them using the provided docker image. More details about how to build the docker image, and about running the benchmark using this image can be found below.

### Running Benchmarks Using Docker

Please remember to install Docker before following the next instructions.

To run our benchmarks using Docker we must first build our image using the Dockerfile provided in this repository. To further simplify this process we also provide the build-img.sh bash script that automatically builds the Docker image. To execute the script simply run ``bash build-img.sh`` while in the same directory as the script ```build-img.sh```. After having built the image, we can execute the image using the ```run-it.sh``` script, which can be executed by running the command ```bash run-it.sh```. We should now find our selves at directory ```/home``` of our container. 

In the ```/home``` directory of our container you should find all of our source code along with build scripts. To build our benchmarks we must simply run the command ```bash build-bench.sh``` while at the ```/home```. This command will build all our benchmark executables, placing them inside ```/home/build```. Here is a list of the Fuzzy PSI protocols and their respective benchmark executable names:

- L Infinity Fuzzy PSI: ```fuzzylinf_bench```
- L1 Fuzzy PSI: ```fuzzyl1_bench```
- L2 Fuzzy PSI: ```fuzzyl2_bench```

All the benchmarks were implemented using the widely known and well documented [Catch2](https://github.com/catchorg/Catch2) C++ library. Each of these executables implement multiple benchmark, each one with different parameters. To list all the benchmarks implemented by an executable, simply run the executable with the parameter ```--list-tests```. For example:

```./fuzzylinf_bench --list-tests```

This will output a list of benchmark names. One example of a benchmark name is ```fuzzyl1 (n=m=256 d=2 delta=10 ssp=40)```. We can run this specific benchmark running the following command:

```./fuzzyl1_bench --benchmark-samples 1 "fuzzyl1 (n=m=256 d=2 delta=10 ssp=40)"```

The extra parameter ```--benchmark-samples``` tells Catch2 how many times to run the benchmark before calculating the final benchmark statistics. In this case we will be running the benchmark just a single time. Here is a Catch2 benchmark output example:

[PUT IMAGE HERE]

