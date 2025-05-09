def read_input(filename="/home/kaua/vscode/ufpi/bioinformatica/SmithWaterman/input.txt"):
    with open(filename, 'r') as file:
        lines = [line.strip() for line in file.readlines()]
    seq1 = lines[0]  # vertical
    seq2 = lines[1]  # horizontal
    gap = int(lines[2])
    mismatch = int(lines[3])
    match = int(lines[4])
    return seq1, seq2, gap, mismatch, match

def score(a, b, match, mismatch):
    return match if a == b else mismatch

def smith_waterman(seq1, seq2, gap, mismatch, match):
    m, n = len(seq1), len(seq2)
    H = [[0 for _ in range(n+1)] for _ in range(m+1)]
    max_score = 0
    max_pos = (0, 0)

    # Fill score matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
            match_score = H[i-1][j-1] + score(seq1[i-1], seq2[j-1], match, mismatch)
            delete = H[i-1][j] + gap
            insert = H[i][j-1] + gap
            H[i][j] = max(0, match_score, delete, insert)

    # encontrar a célula de maior score na última coluna
    max_i = 0
    for i in range(1, m+1):
        if H[i][n] > max_score:
            max_score = H[i][n]
            max_i = i
    max_pos = (max_i, n)

    return H, max_score, max_pos

def backtrace(H, seq1, seq2, gap, mismatch, match, start_pos):
    aligned1, aligned2 = "", ""
    i, j = start_pos
    while i > 0 and j > 0:
        current = H[i][j]
        diagonal = H[i-1][j-1]
        up = H[i-1][j]
        left = H[i][j-1]

        if current == 0:
            break
        elif current == diagonal + score(seq1[i-1], seq2[j-1], match, mismatch):
            aligned1 = seq1[i-1] + aligned1
            aligned2 = seq2[j-1] + aligned2
            i -= 1
            j -= 1
        elif current == up + gap:
            aligned1 = seq1[i-1] + aligned1
            aligned2 = '-' + aligned2
            i -= 1
        elif current == left + gap:
            aligned1 = '-' + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1
    return aligned1, aligned2

def print_matrix(H, seq1, seq2):
    print("-----------------------------------------------------------")
    print("** matrix **")
    print("=" * 60)
    header = "    " + "  ".join(seq2)
    print(header)
    for i in range(len(H)):
        row_label = seq1[i-1] if i > 0 else " "
        row = "  ".join(f"{val}" for val in H[i])
        print(f"{row_label} {row}")
    print("=" * 60)

def main():
    seq1, seq2, gap, mismatch, match = read_input()
    H, max_score, start_pos = smith_waterman(seq1, seq2, gap, mismatch, match)
    aligned1, aligned2 = backtrace(H, seq1, seq2, gap, mismatch, match, start_pos)

    output = []
    output.append("-----------------------------------------------------------")
    output.append("** matrix **")
    output.append("=" * 60)
    header = "    " + "  ".join(seq2)
    output.append(header)
    for i in range(len(H)):
        row_label = seq1[i-1] if i > 0 else " "
        row = "  ".join(f"{val}" for val in H[i])
        output.append(f"{row_label} {row}")
    output.append("=" * 60)
    output.append(f"Score = {max_score}")
    output.append(f"** Match = {match} | mismatch = {mismatch} | Gap = {gap} **")
    output.append("-----------------------------------------------------------")
    output.append("Alinhamento")
    output.append(aligned1)
    output.append(aligned2)
    output.append("-----------------------------------------------------------")

    # Print on screen
    for line in output:
        print(line)

    # Save to output.txt
    with open("output.txt", "w") as f:
        for line in output:
            f.write(line + "\n")

if __name__ == "__main__":
    main()
