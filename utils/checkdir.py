import os


def check_dir(folder):
    if not os.path.exists(folder):
        os.system(f'mkdir {folder}')


