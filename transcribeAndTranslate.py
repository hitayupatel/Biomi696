from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


def transcribe_sequence(sequence):
    """Transcribes the given sequence.
    :param sequence: A string containing sequence
    :return Sequence: Stranded message RNA"""
    coding_dna = Seq(sequence, IUPAC.unambiguous_dna)
    return coding_dna.transcribe()


def translate_rna(message):
    """Translates the given message RNA.
        :param message: A stranded message RNA
        :return protein: returns a protein"""
    return message.translate()


def main():
    sequence1 = "ATGATTGGCCCGGTTTTTTAA"
    sequence2 = "GTGGTGGGGAAATTCCGCTGA"
    print("Sequence 1: ",sequence1)
    transcribed_sequence1 = transcribe_sequence(sequence1)
    translated_sequence1 = translate_rna(transcribed_sequence1)
    print("transcribe_sequence: ", transcribed_sequence1)
    print("translate_sequence: ", translated_sequence1)

    print("Sequence 2: ", sequence2)
    transcribed_sequence2 = transcribe_sequence(sequence2)
    translated_sequence2 = translate_rna(transcribed_sequence2)
    print("transcribe_sequence: ", transcribed_sequence2)
    print("translate_sequence2: ", translated_sequence2)


if __name__ == "__main__":
    main()
