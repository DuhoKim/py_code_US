import glob
from PyPDF2 import PdfFileMerger

pdfs = glob.glob('/Users/dhk/Documents/publish/ngcic_rev2/figset2/figset2_*.pdf')

merger = PdfFileMerger()

for pdf in pdfs:
    merger.append(open(pdf, 'rb'))

with open('/Users/dhk/Documents/publish/ngcic_rev2/figset2/figset2_all.pdf', 'wb') as fout:
    merger.write(fout)