import subprocess
from sys import platform

if platform == "linux":
    str_add = ""
elif platform == "win32" or platform == "win64":
    str_add = ".exe"
else:
    raise OSError("Yeah I dunno what kinda OS you're running but it's wrong :)")
    

#print("Meanline File Generation Complete!")

# Call excecutable stagen-18.1.exe
stagen = subprocess.Popen('Location of stagen-18.1.exe')

print("Stage File Generation Complete!")

# Call excecuable multall-open-20.9.exe
multall = subprocess.Popen('multall-open-20.9.exe')

print("CFD of Turbomachinery Stage Complete!")
