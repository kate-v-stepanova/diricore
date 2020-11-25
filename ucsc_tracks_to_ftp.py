import ftplib
import glob
import os
import datetime
from datetime import datetime
import sys
from pathlib import Path

project = None
if len(sys.argv) >= 2:
    project = sys.argv[1]

bam_type = "hq_unique"
if len(sys.argv) >= 3:
    bam_type = sys.argv[2] 

remove_files = False
if len(sys.argv) >= 4:
    if sys.argv[3] == 'remove':
        remove_files = True

host = "ftp.dkfz-heidelberg.de"
BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
#tracks_dir = "analysis/output/gen_tracks" #.format(bam_type) #ucsc_track_annotation.txt"
ftp = ftplib.FTP()
ftp.connect(host)
ftp.login()

ftppath = '/outgoing/B250/ucsc'
print(project)
now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print(now)

def mkdirs(path):
    try:
        ftp.cwd(path)
    except:
        if path == ftppath:
            print('Path is {}'.format(ftppath))
            return
        try:
            ftp.mkd(path)
        except:
            mkdirs(str(Path(path).parent))
            mkdirs(path)

def remove(project_id, bam_type=""):
    proj_path = os.path.join(ftppath, project_id, bam_type)
    ftp.cwd(proj_path)
    for track in ftp.nlst():
        ftp.delete(track)



def upload_project(project_id, bam_type=""):
   print(project_id, bam_type)
   tracks_dir = os.path.join(BASE_DIR, project_id, "analysis/output/gen_tracks/{}".format(bam_type))
   tracks_pattern = os.path.join(tracks_dir, "*.bw")
   annotation_path = os.path.join(tracks_dir, "ucsc_track_annotation.txt")
   #for project_path in glob.glob(BASE_DIR + '/*/'):
   # project_path will be: /icgc/dkfzlsdf/analysis/OE0532/3784/
   # removing last slash
   #project_id = os.path.basename(project_path[:-1]) 
   project_path = os.path.join(BASE_DIR, project_id)
   if os.path.isfile(os.path.join(annotation_path)):
#   if os.path.isfile(os.path.join(project_path, tracks_dir, bam_type, "ucsc_track_annotation.txt")):
       ftp_proj_path = os.path.join(ftppath, project_id, bam_type)
       try:
           ftp.cwd(ftp_proj_path)
           ftp_files = ftp.nlst(ftp_proj_path)
           if len(ftp_files) != len(glob.glob(tracks_pattern)):
               for ucsc_track in glob.glob(tracks_pattern):
                   f = open(ucsc_track, "rb")
                   ftp.storbinary('STOR ' + os.path.basename(ucsc_track), f)
       except:
           print(os.path.dirname(ftp_proj_path))
           try:
               mkdirs(ftp_proj_path)
#               ftp.mkd(os.path.dirname(ftp_proj_path))
           except:
               print('Cant create path: {}'.format(ftp_proj_path))
           print("Uploading files to: {}".format(ftp_proj_path))
           for ucsc_track in glob.glob(tracks_pattern):
               print(os.path.basename(ucsc_track))
               f = open(ucsc_track, "rb")
               ftp.storbinary('STOR ' + os.path.basename(ucsc_track), f)
       ftp.cwd(ftppath)
   else:
        print('Annotation not found: {}'.format(annotation_path))

if project is None:
    for project_path in glob.glob(BASE_DIR+ "/*/"):
        project = os.path.basename(project_path[:-1])
        print(project)
        for bam_type in ['hq', 'hq_unique', 'all', 'all_unique', ""]:
            upload_project(project, bam_type)
elif remove_files:
    print('Removing tracks')
    remove(project, bam_type)
else:
    print('Uploading tracks')
    upload_project(project, bam_type)

ftp.quit()    
