class FileUtils:
    def __init__(self, name, path):
        self.name = name
        self.path = path

    def read_file(self) -> str:
        """
        Function: Read Gen File
        Purpose: Reads a genomic sequence file (FASTA-like format), removes
        whitespace and ignores metadata lines (those starting with '<'),
        then concatenates the sequence into a single string.

        Complexity: O(n), where n is the number of characters in the file.
        """
        try:
            with open(self.path, "r", encoding="utf-8") as f:
                seq = [
                    line.strip()
                    for line in f
                    if not line.startswith("<")
                ]
            return "".join(seq)
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {self.path}")
        except OSError as e:
            raise RuntimeError(f"Error reading file {self.path}: {e}")

    def read_protein_file(self) -> dict:
        """
        Function: Read Protein File
        Purpose: Reads a protein FASTA file and stores each protein sequence
        in a dictionary, where the key is the protein name (header line without '>'),
        and the value is the corresponding amino acid sequence.

        Complexity: O(n), where n is the number of characters in the file.
        """
        protein_dict = {}
        current_sequence = ""
        try:
            with open(self.path, "r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith(">"):
                        current_name = line[1:].strip()
                    else:
                        current_sequence = line.strip()
                        protein_dict[current_name] = current_sequence
            return protein_dict
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {self.path}")
        except OSError as e:
            raise RuntimeError(f"Error reading file {self.path}: {e}")


class Utils:

    @staticmethod
    def menu() -> str:
        """
        Function: Main menu
        Pourpose:
        Complexity: 
        """
        ...

    @staticmethod
    def list_files():
        archivos = [
            "gen-M.txt",
            "gen-ORF1AB.txt",
            "gen-S.txt",
            "SARS-COV-2-MN908947.3.txt",
            "SARS-COV-2-MT106054.1.txt",
            "seq-proteins.txt"
        ]
        print("Archivos disponibles:")
        for i in archivos:
            print(f"---{archivos[i]}---")
        

    @staticmethod
    def format_manchester(s: str) -> str:
        """
        Function: Format list
        Pourpose: Give format to the list in order to be processed by
        Manchester alglrithm in algorihms.py file.
        Complexity: O(n)
        """
        simbol_1 = '@'
        simbol_2 = '#'
        simbol_3 = '$'
        char = 0

        result = [simbol_1]
        new_size = 2 * len(s) + 3

        for i in range(1, new_size - 1):
            if i % 2 != 0:
                result.append(simbol_3)
            else:
                if char < len(s):
                    result.append(s[char])
                    char += 1

        result.append(simbol_2)
        return result