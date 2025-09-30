from utils import Utils

class Algorithms:

    @staticmethod
    def mancher_algorithm(s: str) -> str:
        """
        Function: Manacher's Algorithm
        Purpose: Find the largest palindrome substring in a given string.
        Complexity: O(n), where n is the length of the input string.
        """
        simbol_1 = '@'
        simbol_2 = '#'
        simbol_3 = '$'

        new_list = Utils.format_manchester(s)  # formateador
        n = len(new_list)
        if n == 0:
            return ""

        P = [0] * n
        center = 0
        limit = 0

        for i in range(1, n - 1):
            if i < limit:
                simetric = 2 * center - i
                P[i] = min(limit - i, P[simetric])

            gap = P[i] + 1
            # expandir cuidando los límites
            while i - gap >= 0 and i + gap < n and new_list[i - gap] == new_list[i + gap]:
                P[i] += 1
                gap += 1

            if i + P[i] > limit:
                limit = i + P[i]
                center = i

        # obtener el mejor centro
        best_center = max(range(n), key=lambda k: P[k])
        best_radius = P[best_center]
        L = best_center - best_radius
        R = best_center + best_radius

        # reconstruir
        out = []
        for k in range(L, R + 1):
            ch = new_list[k]
            if ch not in (simbol_1, simbol_2, simbol_3):
                out.append(ch)

        return "".join(out)

    @staticmethod
    def computeLPSArray(pat, M, lps):
        """
        Function: Compute LPS (Longest Prefix Suffix) Array
        Purpose: Preprocess the pattern to generate the LPS array used by the 
                 Knuth-Morris-Pratt (KMP) algorithm. LPS array stores the length 
                 of the longest prefix which is also a suffix for each sub-pattern.
        Complexity: O(M), where M is the length of the pattern.
        """
        len = 0
        i = 1
        lps[0] = 0
        while i < M:
            if pat[i] == pat[len]:
                lps[i] = len + 1
                len += 1
                i += 1
            else:
                if len != 0:
                    len = lps[len - 1]
                else:
                    lps[i] = 0
                    i += 1 

    @staticmethod
    def KMPSearch(pat, txt, nombre):
        """
        Function: KMP Search Algorithm
        Purpose: Search for all occurrences of a pattern (pat) inside a text (txt),
                 using the preprocessed LPS array to avoid redundant comparisons.
                 Prints the starting index of each match along with the first 12
                 nucleotides of the match.
        Complexity: O(N + M), where N is the length of the text and M is the length of the pattern.
        """
        N = len(txt)
        M = len(pat)
        lps = [0]*M
        Utils.computeLPSArray(pat, M, lps)
        i = 0
        j = 0
        while i < N:
            if txt[i] == pat[j]:
                i += 1
                j += 1
            if j == M:
                print("El inicio del genoma ", nombre," es", i - j)
                print(txt[i-j:i-j + 12]) #Se imprimen los primeros 12 nucleotidos del patron encontrado
                j = lps[j-1]
            elif i < N and txt[i] != pat[j]:
                if j != 0:
                    j = lps[j-1]
                else:
                    i += 1
