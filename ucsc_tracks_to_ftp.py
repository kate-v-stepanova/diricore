import ftplib
import glob
import os
import datetime
from datetime import datetime
import sys

project = None
if len(sys.argv) == 2:
    project = sys.argv[1]

host = "ftp.dkfz-heidelberg.de"
BASE_DIR = "/icgc/dkfzlsdf/analysis/OE0532"
tracks_dir = "analysis/output/gen_tracks" #ucsc_track_annotation.txt"
ftp = ftplib.FTP()
ftp.connect(host)
ftp.login()

file = '/home/rock/test.txt'
ftppath = '/outgoing/B250/ucsc'

now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print(now)

def upload_project(project_id):
  #for project_path in glob.glob(BASE_DIR + '/*/'):
   # project_path will be: /icgc/dkfzlsdf/analysis/OE0532/3784/
   # removing last slash
   #project_id = os.path.basename(project_path[:-1]) 
   project_path = os.path.join(BASE_DIR, project_id)
   if os.path.isfile(os.path.join(project_path, tracks_dir, "ucsc_track_annotation.txt")):
       proj_path = os.path.join(ftppath, project_id)
       try:
           ftp.cwd(proj_path)
           files = ftp.nlst(proj_path)
           if len(files) == 0:
               for ucsc_track in glob.glob(os.path.join(project_path, tracks_dir, "*.bw")):
                   f = open(ucsc_track, "rb")
                   ftp.storbinary('STOR ' + os.path.basename(ucsc_track), f)
       except:
           ftp.mkd(proj_path)
           ftp.cwd(proj_path)
           print("Uploading files to: {}".format(proj_path))
           for ucsc_track in glob.glob(os.path.join(project_path, tracks_dir, "*.bw")):
               print(os.path.basename(ucsc_track))
               f = open(ucsc_track, "rb")
               ftp.storbinary('STOR ' + os.path.basename(ucsc_track), f)
       ftp.cwd(ftppath)
   else:
        for ucsc_dir in glob.glob(os.path.join(project_path, tracks_dir,'*') + os.path.sep): # dirs only
            print(ucsc_dir)
            ucsc_annotation = os.path.join(ucsc_dir, "ucsc_track_annotation.txt")
            bam_type = os.path.basename(ucsc_dir[:-1]) # remove last /
            proj_path = os.path.join(ftppath, project_id)
            try:
                ftp.cwd(proj_path)
            except:
                ftp.mkd(proj_path)
                proj_path = os.path.join(proj_path, bam_type)
                ftp.mkd(proj_path)
                ftp.cwd(proj_path)
                print("Uploading files to: {}".format(proj_path))
                for ucsc_track in glob.glob(os.path.join(ucsc_dir, "*.bw")):
                    print(os.path.basename(ucsc_track))
                    f = open(ucsc_track, "rb")
                    ftp.storbinary("STOR " + os.path.basename(ucsc_track), f)
                ftp.cwd(ftppath)

if project is None:
    for project_path in glob.glob(BASE_DIR+ "/*/"):
        project = os.path.basename(project_path[:-1])
        print(project)
        upload_project(project)
else:
    upload_project(project)

ftp.quit()    
