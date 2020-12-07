import re


def filereader(filename):
    """Leest een gegeven file in ??? format in. Slaat de markers met
    bijbehorende allelen op in een dictionary. Markers die niet door de
    chi-squared test komen worden niet opgeslagen.

    Input:  filename: str, path naar input bestand.

    Output: marker_dict: dict, dictionary met markers en allelen.
            total: int, aantal markers in marker_dict.
    """
    marker_dict = {}
    marker = ""
    alleles = []
    with open(filename) as file:
        # Overslaan eerste 7 regels.
        for i in range(7):
            line = file.readline()

        for line in file:
            if not line.startswith(" "):
                if marker:
                    # print(marker, chi_squared(alleles))
                    # Slaat markers alleen op als ze door chi_squared komen.
                    if chi_squared(alleles, marker):
                        marker_dict[marker] = alleles
                    # Reset lijst met allelen.
                    alleles = []
                marker = line.split()[0]
            else:
                # Splitst regel in een lijst met allelen.
                alleles.extend(line.split())
    if chi_squared(alleles, marker):
        marker_dict[marker] = alleles
    total = len(marker_dict[marker])

    return marker_dict, total


def compare(marker_dict):
    """Compares all markers with eachother and saves the amount of
    non-matches alleles to a dictionary.

    Input:  marker_dict - dict, dictionary containing markers and alleles.

    Output: score_dict - dict, dictionary containing every comparison
                            and the amount of mismatches.
    """

    score_dict = {}
    # Iterates over all markers.
    for i in range(len(marker_dict)):
        marker1 = list(marker_dict.keys())[i]
        marker1_list = marker_dict[marker1]
        # Iterates over all following markers.
        for j in range(i + 1, len(marker_dict)):
            marker2 = list(marker_dict.keys())[j]
            marker2_list = marker_dict[marker2]
            som = 0
            # Iterates over the alleles simutaneously and records the
            # amount of mismatches (excluding missing values).
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
    """Returns combination frequency."""
    return (score / total) * 100


def get_factors(score_dict, total):
    """Returns converts scores in a dictionary to recombination frequencies."""
    for vergelijking, score in score_dict.items():
        score_dict[vergelijking] = calc_recombinance(score, total)
    return score_dict


def chi_squared(data: list, marker: str, output_file="Chi-output.txt"):
    """Voert chi-squared test uit op een gegeven lijst met allelen 'a' & 'b'.
    (df = 1, a = 0.05). Uitkomst wordt ook toegevoegd aan output_file.

    Input:  data - list, list containing alleles..
            marker - str, marker name.
            output_file - str, path for the output file.
                (default: 'Chi-output.txt')

    Output: True if marker is usable, otherwise False.
    """
    count_a = data.count("a")
    count_b = data.count("b")
    expected_a = len(data) / 2
    expected_b = len(data) / 2

    temp_a = ((count_a - expected_a)**2) / expected_a
    temp_b = ((count_b - expected_b)**2) / expected_b

    outcome = temp_a + temp_b

    with open(output_file, "a+") as file:
        file.write(f"{marker}, a: {count_a}, b: {count_b}, "
                   f"expected: {expected_a}, total: {len(data)}, "
                   f"outcome: {outcome}, accepted: {outcome <= 3.84}\n")

    # Returnt True als outcome <= 3.84 (df = 1, a = 0.05)
    if outcome <= 3.84:
        return True
    else:
        return False


def get_afstanden(factors):
    """Looks for the two markers with the highest recombination
    frequency. Chooses one as a starting point, then searches al
    comparisons with that marker and sorts them based on frequency.
    """
    # Looks for the markers with the highest recombination frequency.
    start_marker = sorted(tuple(factors.items()),
                          key=lambda x: x[1], reverse=True)[0][0].split(",")[0]

    matches = []
    # Looks for all recombination frequncies between the start_marker
    # and other markers.
    for key, value in factors.items():
        if re.search(start_marker, key):
            temp_list = key.split(", ")
            if temp_list[0] != start_marker:
                matches.append((temp_list[0], value))
            else:
                matches.append((temp_list[1], value))

    # Adds start_marker.
    matches.append((start_marker, 0))

    # Sorts markers by distance.
    matches.sort(key=lambda x: x[1])
    print(len(matches))
    return matches


def write_file(matches):
    """Writes the distances between markers to a file in .tsv format."""
    with open("output.csv", "w") as file:
        file.write("Group\t1\n")

        for match in matches:
            file.write(f"{match[0]}\t{match[1]}\n")


def main():
    filename = "CvixLer-MarkerSubset-LG1.txt"
    marker_dict, total = filereader(filename)
    score_dict = compare(marker_dict)
    print(score_dict)
    factors = get_factors(score_dict, total)
    matches = get_afstanden(factors)
    write_file(matches)
    print()


main()
