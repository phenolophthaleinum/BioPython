from Bio.SeqIO import parse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna
from Bio import Alphabet, Entrez
from Bio.SeqUtils import GC
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import motifs

__author__ = 'Maciej Michalczyk'
__version__ = '09122019'

class CurrentSession:

    def __init__(self, sequence = None, comp_seq = None, transcribed_seq = None,
                 translated_seq = None, id = None, name = None, desc = None, gc_perc = None,
                 record = None):
        self.sequence = None
        self.comp_seq = None
        self.transcribed_seq = None
        self.translated_seq = None
        self.id = None
        self.name = None
        self.desc = None
        self.gc_perc = None
        self.record = SeqRecord(self.sequence, id = self.id)
        self.file_session = None


    def openFile(self, filename):
        file = open(filename + '.fasta')
        self.file_session = str(filename + '.fasta')
        return file


    def closeFile(self, file_handle):
        file_handle.close()


    def getSequenceInfo(self, file_handle):
        records = parse(file_handle, "fasta")

        for record in records:
            self.id = record.id
            self.name = record.name
            self.desc = record.description
            self.sequence = Seq(str(record.seq), IUPAC.ambiguous_dna)
            print("ID: {}".format(self.id))
            print("Name: {}".format(self.name))
            print("Description: {}".format(self.desc))
            print("Sequence: {}".format(self.sequence))
            # print("Complementary sequence: {}".format(sequence.complement()))
            print("------------------------------------------------------------")
        return


    def getComplementarySequence(self, file_handle):
        records = parse(file_handle, "fasta")

        for record in records:
            self.sequence = Seq(str(record.seq), IUPAC.unambiguous_dna)
            self.name = record.name
            self.comp_seq = self.sequence.complement()
            print("Name: {}".format(self.name))
            if Alphabet._verify_alphabet(self.sequence) == True:
                print("Sequence: {}".format(self.sequence))
                print("Complementary sequence: {}".format(self.comp_seq))
                print("------------------------------------------------------------")
            else:
                print("This sequence is not a DNA, can't get a complementary of that. Load correct sequence.")
        return


    def transcribeSequence(self, file_handle):
        records = parse(file_handle, "fasta")

        for record in records:
            self.sequence = Seq(str(record.seq), IUPAC.unambiguous_dna)
            self.name = record.name
            self.transcribed_seq = self.sequence.transcribe()
            print("Name: {}".format(self.name))
            if Alphabet._verify_alphabet(self.sequence) == True:
                print("Sequence: {}".format(self.sequence))
                print("Transcribed sequence: {}".format(self.transcribed_seq))
                print("------------------------------------------------------------")
            else:
                print("This sequence is not a DNA, can't get a complementary of that. Load correct sequence.")


    def translateSequence(self, file_handle, stop):
        records = parse(file_handle, "fasta")

        for record in records:
            self.sequence = Seq(str(record.seq), IUPAC.unambiguous_rna)
            self.name = record.name
            print("Name: {}".format(self.name))
            if Alphabet._verify_alphabet(self.sequence) == True and stop == 'y':
                self.translated_seq = self.sequence.translate(to_stop=True)
                print("Sequence: {}".format(self.sequence))
                print("Translated sequence: {}".format(self.translated_seq))
                print("------------------------------------------------------------")
            elif Alphabet._verify_alphabet(self.sequence) == True and stop == 'n':
                self.translated_seq = self.sequence.translate()
                print("Sequence: {}".format(self.sequence))
                print("Translated sequence: {}".format(self.translated_seq))
                print("------------------------------------------------------------")
            else:
                print("This sequence is not a RNA, can't translate that. Load correct sequence.")


    def get_GC_Content(self, file_handle):
        records = parse(file_handle, "fasta")

        for record in records:
            self.sequence = Seq(str(record.seq), IUPAC.unambiguous_dna)
            self.name = record.name
            self.gc_perc = GC(self.sequence)
            print("Name {}".format(self.name))
            if Alphabet._verify_alphabet(self.sequence) == True:
                print("Sequence: {}".format(self.sequence))
                print("GC content: {}%".format(self.gc_perc))
                print("------------------------------------------------------------")
            else:
                print("This sequence is not a DNA, only calculate GC content in DNA. Load correct sequence.")


    def fetchRecord(self, db, accession):
        Entrez.email = "A.N.Other@example.com"
        handle = Entrez.efetch(db = db, id = accession, rettype = "fasta")
        #print(handle.read())

        record = SeqIO.read(handle, "fasta")
        print(record)
        filename = record.name
        file = open(filename + ".fasta", "w")
        SeqIO.write(record, file, "fasta")
        return filename


    def runBlast(self, type, database):
        seq_record = next(SeqIO.parse(open(self.file_session), 'fasta'))
        print("Requesting BLAST (might take a few minutes...)")
        request_handle = NCBIWWW.qblast(type, database, seq_record.seq)
        print("BLAST succeeded.")

        with open("{}_blast.xml".format(self.file_session), "w") as save_f:
            save_f.write(request_handle.read())
        request_handle.close()
        print("BLAST results saved.")


    def alignPairwise(self, file_handle, alignment_type):
        try:
            records = parse(file_handle, "fasta")
            number = 1
            seq1 = None
            seq2 = None

            for record in records:
                if number == 1:
                    seq1 = record.seq
                elif number == 2:
                    seq2 = record.seq
                number += 1
            if seq2 is None:
                print("Error: There is only one sequence in the file.")
                return

            if alignment_type == str(1):
                alignments = pairwise2.align.globalxx(seq1, seq2)
            elif alignment_type == str(2):
                alignments = pairwise2.align.localxx(seq1, seq2)
            elif alignment_type == str(3):
                match = int(input("Define points given for match: "))
                mismatch = int(input("Define points deduced for mismatch: "))
                o_gap = int(input("Define penalty for gap opening: "))
                ext_gap = int(input("Define penalty for gap extension: "))
                alignments = pairwise2.align.globalms(seq1, seq2, match, mismatch,
                                                          o_gap, ext_gap)

            for alignment in alignments:
                print("RAW ALIGNMENT: ")
                print(alignment)
                print("FORMATTED ALIGNMENT: ")
                print(format_alignment(*alignment))
        except Exception as e1:
            print("Error, problably there is only one sequence.")


    def createMotif(self, file_handle):
        records = parse(file_handle, "fasta")
        logofile = self.file_session + "_logo.png"
        seqs_motif = []

        for record in records:
            self.sequence = Seq(str(record.seq))
            seqs_motif.append(self.sequence)
        seqs = motifs.create(seqs_motif)
        print(seqs.counts)
        seqs.weblogo(logofile)
        print("Weblogo saved.")

    def getElems(self, length):
        for base in range(len(self.sequence)):
            if base + length > len(self.sequence):
                break
            else:
                yield self.sequence[base:base + length]

    def saveActions(self):
        if self.file_session != None:
            with open("{}_actions.txt".format(self.id), "w") as save_f:
                save_f.write("""ID: {}
NAME: {}
DESCRIPTION: {}
ORIGINAL SEQUENCE: {}
COMPLEMENTARY SEQUENCE: {}
TRANSCRIBED SEQUENCE: {}
TRANSLATED SEQUENCE: {}
G/C PERCENTAGE: {}%""".format(self.id, self.name, self.desc, self.sequence, self.comp_seq,
                              self.transcribed_seq, self.translated_seq, self.gc_perc))
            save_f.close()
            print("Your actions were saved!")
        else:
            print("Nothing to save, probably you haven't loaded any file before.")
   # def convertFASTAtoGENBANK(self, filename):
   #     file = open(filename + ".fasta")

   #     record = SeqIO.read(file, "fasta")
   #     record.seq.alphabet = generic_dna

   #     file_genbank = open(filename + ".gbk", "w")
   #     SeqIO.write(record, file_genbank, "genbank")
   #     file_genbank.close()
   #     file.close()
if __name__ == '__main__':

    session = CurrentSession()

    while True:
        print("""////Bio Python////
1. Load FASTA file
2. Load record info
    3. Get complementary sequence
    4. Transcribe sequence
    5. Translate sequence
    6. Get GC content
7. Fetch and load FASTA from Entrez
*8. Convert FASTA to GenBank
    9. Run BLAST
    10. Perform pairwise aligment
    11. Create motifs and weblogo
    12. Save your actions made on FASTA file to txt file
    13. Print sequence substrings
    
=== Current session file: {} ===

Type 'help' for help.
Type 'quit' to exit.""".format(session.file_session))

        menu_pos = input('>>').lower()

        if menu_pos == str(1):
            try:
                print("Type name of FASTA file to process: ")
                filename = input()
                file_handle = session.openFile(filename)
                print("FASTA loaded!")
            except Exception as e:
                print("No such file or directory.")
        elif menu_pos == str(2):
            try:
                file_handle = session.openFile(filename)
                session.getSequenceInfo(file_handle)
                session.closeFile(file_handle)
            except Exception as e:
                print("File is not loaded.")
        elif menu_pos == str(3):
            try:
                file_handle = session.openFile(filename)
                session.getComplementarySequence(file_handle)
                session.closeFile(file_handle)
            except Exception as e:
                print("File is not loaded.")
        elif menu_pos == str(4):
            try:
                file_handle = session.openFile(filename)
                session.transcribeSequence(file_handle)
                session.closeFile(file_handle)
            except Exception as e:
                print("File is not loaded.")
        elif menu_pos == str(5):
            stop = input('Stop translating at first stop codon? [y/n]').lower()
            try:
                file_handle = session.openFile(filename)
                session.translateSequence(file_handle, stop)
                session.closeFile(file_handle)
            except Exception as e:
                print("File is not loaded.")
        elif menu_pos == str(6):
            try:
                file_handle = session.openFile(filename)
                session.get_GC_Content(file_handle)
                session.closeFile(file_handle)
            except Exception as e:
                print("File is not loaded.")
        elif menu_pos == str(7):
            try:
                db = input("Type database name: ").lower()
                accession = input("Type accession to find: ")
                filename = session.fetchRecord(db, accession)
                file_handle = session.openFile(filename)
            except Exception as e:
                print("File is not loaded.")
        elif menu_pos == str(8):
            try:
                print("Type name of FASTA file to process: ")
                filename = input()
               # session.convertFASTAtoGENBANK(filename)
            except Exception as e:
                print("File is not loaded.")
        elif menu_pos == str(9):
            try:
                file_handle = session.openFile(filename)
                type = input("Type the type of BLAST: ")
                database = input("Type database name: ")
                session.runBlast(type, database)
            except Exception as e:
                print("File is not loaded.")
        elif menu_pos == str(10):
            try:
                print("""Choose type of aligment: 
                1. Global Pairwise (default parameters)
                2. Local Pairwise (default parameters)
                3. Global Pairwise with custom parameters""")
                alignment_type = input('>>')
                file_handle = session.openFile(filename)
                session.alignPairwise(file_handle, alignment_type)
                session.closeFile(file_handle)
            except Exception as e:
                print("File is not loaded.")
        elif menu_pos == str(11):
            try:
                file_handle = session.openFile(filename)
                session.createMotif(file_handle)
                session.closeFile(file_handle)
            except Exception as e:
                print("File is not loaded")
        elif menu_pos == str(12):
            session.saveActions()
        elif menu_pos == str(13):
            try:
                length = int(input("Length of substrings:"))
                iterator = session.getElems(length)
                print(session.sequence)
                i = 0
                for base in iterator:
                    print(' ' * i + base)
                    i += 1
                    print(' ' * i + next(iterator))
                    i += 1
            except StopIteration:
                pass
            except Exception as e:
                print("File is not loaded")
        elif menu_pos == 'debug':
            print("{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n".format(session.id, session.name,
                                                            session.desc,session.sequence, 
                                                            session.comp_seq, session.transcribed_seq,
                                                            session.translated_seq, session.gc_perc))
        elif menu_pos == 'quit':
            break
        elif menu_pos == 'help':
            print("""
quickHELP:
Indent operations in menu needs file to be opened first.
Be patient while doing BLAST.
If in menu something is marked with an asterisk, then it is not usable.
Have fun!
""")
        else:
            print("Unknown command.")
