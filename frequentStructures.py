#!  /usr/bin/env python3
import RNA
import random

def generateSeq():
    seq = ""

    for i in range(15):
        seq += random.choice(["A","U","C","G"])

    return seq


def main():
    structures = dict()
    out = []


    for i in range(int(1e7)):
        if i % 1000 == 0:
            print(i)
        seq = generateSeq()

        (ss, mfe) = RNA.fold(seq)

        if ss in structures:
            structures[ss]["count"] += 1

        else:
            structures[ss] = {"count": 1, "sequence": seq }

    for key in structures:
        out.append([key, structures[key]["count"], structures[key]["sequence"]])

    out = sorted(out, key = lambda l: l[1], reverse = True)

    with open("Frequent structure.txt", "w") as f:
        f.write(f'Structure\t\t\t\tCount\tSequence\n')

        for index ,row in enumerate(out):
            if index < 100:
                f.write(f'{row[0]}\t{row[1]}\t\t{row[2]}\n')
            else:
                break
       

if __name__ == "__main__":
    main()
