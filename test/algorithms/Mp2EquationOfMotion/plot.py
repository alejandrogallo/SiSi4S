import numpy
import matplotlib.pyplot as plt
import os

folder = os.path.dirname(__file__)
print(folder)

file_name = folder+"/"+"SimlarityTransformedHamiltonianSD.dat"
data = numpy.loadtxt(file_name, skiprows=2)
# print(data)
plt.matshow(data)
fig = plt.gcf()
plt.colorbar()


# print("Determinant")
try:
    det = numpy.linalg.det(data)
except:
    det = 0

if det:
    eigva, eigve = numpy.linalg.eig(data)

    print("\n\nEigenvalues")
    for e in eigva:
        print(e)

    # print("\n\nEigenvalues(sorted)")
    # for e in eigva.sort():
        # print(e)

    print("\n\nEigenvectors")
    for ev in eigve:
        print(ev)

plt.savefig(folder+"/matrix.pdf")
plt.savefig("matrix.pdf")

#vim-run: python % && mupdf matrix.pdf
#vim-run: python %
