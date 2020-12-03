def filereader(filename):
    marker_dict = {}
    marker = ""
    alleles = []
    with open(filename) as file:
        for i in range(7):
            line = file.readline()

        for line in file:
            if not line.startswith(" "):
                if marker:
                    marker_dict[marker] = alleles
                    alleles = []
                marker = line.split()[0]
            else:
                alleles.extend(line.split())
    marker_dict[marker] = alleles

    print(marker_dict)


def main():
    filename = "CvixLer-MarkerSubset-LG1.txt"
    filereader(filename)



main()
