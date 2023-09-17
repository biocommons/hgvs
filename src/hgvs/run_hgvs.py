import subprocess
# subprocess: used to interact with the system's shell or command-line environment.
class HGVSShellLauncher:
    def __init__(self):
        self.command = 'hgvs-shell'

    def start_shell(self):
        try:
            subprocess.run(self.command, shell=True)
        except Exception as e:
            print(f"An error occurred: {e}")

if __name__ == "__main__":
    shell_launcher = HGVSShellLauncher()
    shell_launcher.start_shell()
