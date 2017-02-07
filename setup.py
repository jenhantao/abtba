import os
# download mouse genome from UCSC
os.system('wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz')
# decomporess genome
os.system('tar -zxvf chromFA.tar.gz')
os.mkdir('./mm10')
os.system('mv ./*fa ./mm10')
