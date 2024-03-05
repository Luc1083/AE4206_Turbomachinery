import subprocess

# Call excecutable meangen-17.4.exe
meangen = subprocess.popen('Location of meangen-17.4.exe')
meangen.stdin.write('F')

print("Meanline File Generation Complete!")

# Call excecutable stagen-18.1.exe
stagen = subprocess.popen('Location of stagen-18.1.exe')

print("Stage File Generation Complete!")

# Call excecuable multall-open-20.9.exe
multall = subprocess.popen('multall-open-20.9.exe')

print("CFD of Turbomachinery Stage Complete!")
