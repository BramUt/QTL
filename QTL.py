import re

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
    total = len(marker_dict[marker])

    # print(marker_dict)
    return marker_dict, total


def compare(marker_dict):
    score_dict = {}
    for i in range(len(marker_dict)):
        marker1 = list(marker_dict.keys())[i]
        marker1_list = marker_dict[marker1]
        for j in range(i + 1, len(marker_dict)):
            marker2 = list(marker_dict.keys())[j]
            marker2_list = marker_dict[marker2]
            som = 0
            for k in range(len(marker1_list)):
                allel1 = marker1_list[k]
                allel2 = marker2_list[k]
                if allel1 == "-" or allel2 == "-":
                    continue
                elif allel1 != allel2:
                    som += 1
            score_dict[f"{marker1}, {marker2}"] = som
    return score_dict


def calc_recombinance(score, total):
    return (score / total) * 100

def get_factors(score_dict, total):
    for vergelijking, score in score_dict.items():
        score_dict[vergelijking] = calc_recombinance(score, total)
    return score_dict

def get_afstanden(factors):
    start_marker = sorted(tuple(factors.items()),
                          key=lambda x: x[1], reverse=True)[0][0].split(",")[0]
    # print(sorted(tuple(factors.items()),
    #              key=lambda x: x[1]))
    # print(factors)
    # print(tuple(factors.keys()))
    matches = []
    for key, value in factors.items():
        if re.search(start_marker, key):
            temp_list = key.split(", ")
            if temp_list[0] != start_marker:
                matches.append((temp_list[0], value))
            else:
                matches.append((temp_list[1], value))

    matches.append((start_marker, 0))

    matches.sort(key=lambda x: x[1])
    print(len(matches))
    return matches


def write_file(matches):
    with open("output.csv", "w") as file:
        file.write("Group,1\n")

        for match in matches:
            file.write(f"{match[0]},{match[1]}\n")


def main():
    filename = "CvixLer-MarkerSubset-LG1.txt"
    marker_dict, total = filereader(filename)
    score_dict = compare(marker_dict)
    factors = get_factors(score_dict, total)
    matches = get_afstanden(factors)
    write_file(matches)
    print()



main()
