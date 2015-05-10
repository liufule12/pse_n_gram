__author__ = 'Fule Liu'

import subprocess


if __name__ == "__main__":
    subprocess.Popen("gist-train-svm.exe -train sample.mtx.txt -class sample.labels.txt > sample.weights",
                     shell=True).wait()
    subprocess.Popen("gist-classify -train sample.mtx.txt -learned sample.weights -test test.mtx.txt > test.predict",
                     shell=True).wait()
    print("End")
