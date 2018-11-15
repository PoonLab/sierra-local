# Installing sierra-local

## For Linux Users
The recommended method to install *sierra-local* on a Linux system is to clone the source from this repository and run the `setup.py` installation script with superuser privileges, which will install the package to `/usr/local/lib`:
```console
art@Jesry:~/git$ git clone http://github.com/PoonLab/sierra-local
Cloning into 'sierra-local'...
warning: redirecting to https://github.com/PoonLab/sierra-local/
remote: Enumerating objects: 71, done.
remote: Counting objects: 100% (71/71), done.
remote: Compressing objects: 100% (52/52), done.
remote: Total 709 (delta 40), reused 39 (delta 19), pack-reused 638
Receiving objects: 100% (709/709), 10.08 MiB | 5.54 MiB/s, done.
Resolving deltas: 100% (403/403), done. 
art@Jesry:~/git$ cd sierra-local
art@Jesry:~/git/sierra-local$ sudo python3 setup.py install
[sudo] password for art: 
running install
running build
running build_py
running build_scripts
running install_lib
creating /usr/local/lib/python3.6/dist-packages/sierralocal
copying build/lib/sierralocal/utils.py -> /usr/local/lib/python3.6/dist-packages/sierralocal
[...]
running install_scripts
copying build/scripts-3.6/sierralocal -> /usr/local/bin
changing mode of /usr/local/bin/sierralocal to 775
Changing permissions of /usr/local/lib/python3.6/dist-packages/sierralocal/data to 0o777
art@Jesry:~/git/sierra-local$ sierralocal
usage: sierralocal [-h] [-o OUTFILE] [-xml XML] [--cleanup] [--forceupdate]
                   fasta [fasta ...]
sierralocal: error: the following arguments are required: fasta
```
However, we appreciate that many users will not have root privileges on their computer or may prefer not to clone the git repository, so we provide some alternative installation instructions below.

### Using pip3
To install *sierra-local* from the [Python Package Index](https://pypi.org/), you need to have `pip3` on your system.  You can check whether this program is installed by invoking it from the command line:
```console
art@Jesry:~$ pip3

Command 'pip3' not found, but can be installed with:

sudo apt install python3-pip
```
Note that this response will vary with Linux distributions.  This example uses Ubuntu 18.04.

Alternatively, you can use the `which` command:
```console
art@Jesry:~$ which pip3
art@Jesry:~$
```
If this command returns nothing, then you need to install `pip3` with your system's package manager.  

For Ubuntu/Debian systems, for example:
```console
art@Jesry:~$ sudo apt install python3-pip
[sudo] password for art: 
Reading package lists... Done
Building dependency tree       
Reading state information... Done
The following additional packages will be installed:
  dh-python libexpat1-dev libpython3-dev libpython3.6 libpython3.6-dev libpython3.6-minimal
  libpython3.6-stdlib python-pip-whl python3-dev python3-setuptools python3-wheel python3.6
  python3.6-dev python3.6-minimal
Suggested packages:
  python-setuptools-doc python3.6-venv python3.6-doc binfmt-support
The following NEW packages will be installed:
  dh-python libexpat1-dev libpython3-dev libpython3.6-dev python-pip-whl python3-dev python3-pip
  python3-setuptools python3-wheel python3.6-dev
The following packages will be upgraded:
  libpython3.6 libpython3.6-minimal libpython3.6-stdlib python3.6 python3.6-minimal
5 upgraded, 10 newly installed, 0 to remove and 189 not upgraded.
Need to get 53.1 MB of archives.
After this operation, 82.2 MB of additional disk space will be used.
Do you want to continue? [Y/n] y
```
If you are okay with these packages being installed, then respond with `y`.  Your package manager will retrieve these packages over your network connection and install them in the required locations of your filesystem.  

Now you have two options for installing *sierra-local* with `pip3`:
1. You can install it locally in your user home directory:
   ```console
   art@Jesry:~$ pip3 install sierralocal
   Collecting sierralocal
     Downloading https://files.pythonhosted.org/packages/ce/a8/2501b1b3ad9b8abc7994b7ece3d325b3d765401c4ce7f5760bda765f5267/sierralocal-0.0.1.tar.gz (1.7MB)
       100% |████████████████████████████████| 1.7MB 982kB/s 
   Building wheels for collected packages: sierralocal
     Running setup.py bdist_wheel for sierralocal ... done
     Stored in directory: /home/art/.cache/pip/wheels/e1/2f/28/fab709a2a25b279a11774fbae50520a3dfeb370ad940e84c68
   Successfully built sierralocal
   Installing collected packages: sierralocal
   Successfully installed sierralocal-0.0.1
   ```
   Note that the version number (*e.g., `0.0.1`) will vary over time.
   However, you may find that calling the main program `sierralocal` fails to run anything:
   ```console
   art@Jesry:~$ sierralocal
   sierralocal: command not found
   ```
   This happens because the package was installed in a hidden directory within your home directory:
   ```console
   art@Jesry:~$ locate sierralocal
   /home/art/.local/bin/sierralocal
   ```
   The `.local/` directory is probably not on your `PATH` (the set of directories where the OS looks for executable programs):
   ```console
   art@Jesry:~$ echo $PATH
   /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin
   ```
   To append this directory to the `PATH` environmental variable, you need to run this command:
   ```
   art@Jesry:~$ export PATH=/home/art/.local/bin:$PATH
   art@Jesry:~$ echo $PATH
   /home/art/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin
   ```
   Note this is a temporary fix, and your `PATH` may revert to the original setting.
   
   Now your OS will be able to find the `sierralocal` executable file when you invoke it on the command line:
   ```console
   art@Jesry:~$ sierralocal RT.fa
   searching path /home/art/.local/lib/python3.6/site-packages/sierralocal/data/HIVDB*.xml
   Error: could not find local copy of HIVDB XML, attempting download...
   Updated HIVDB XML from https://hivdb.stanford.edu/assets/media/HIVDB_8.7.9e470b87.xml into /home/art/.local/lib/python3.6/site-packages/sierralocal/data/HIVDB_8.7.9e470b87.xml
   /home/art/.local/lib/python3.6/site-packages/sierralocal/data/apobec.tsv
   Error: could not retrieve APOBEC DRM data
   Updated APOBEC DRMs from https://hivdb.stanford.edu/assets/media/apobec-drms.5b7e1215.tsv into /home/art/.local/lib/python3.6/site-packages/sierralocal/data/apobec.tsv
   HIVdb version 8.7
   Found NucAmino binary /home/art/.local/lib/python3.6/site-packages/sierralocal/bin/nucamino-linux-amd64
   Aligned RT.fa
   100 sequences found in file RT.fa.
   Writing JSON to file RT_results.json
   Time elapsed: 7.9871 seconds (17.795 it/s)
   ```
   Note the `RT.fa` file in this example is a small random sample of HIV-1 RT sequences that was generated using the `retrieve_hivdb_data.py` provided in the `scripts/` directory of this repository.
   
2. If you **do** have root privileges on your computer and you still prefer to use `pip3`, you can do so:
   ```console
   art@Jesry:~$ sudo pip3 install sierralocal
   [sudo] password for art: 
   The directory '/home/art/.cache/pip/http' or its parent directory is not owned by the current user and the cache has been disabled. Please check the permissions and owner of that directory. If executing pip with sudo, you may want sudo's -H flag.
   The directory '/home/art/.cache/pip' or its parent directory is not owned by the current user and caching wheels has been disabled. check the permissions and owner of that directory. If executing pip with sudo, you may want sudo's -H flag.
   Collecting sierralocal
     Downloading https://test-files.pythonhosted.org/packages/1a/67/4b5c97b7bd526bb1b654f9d5fa913f51e70b6c3edb1d2aaefb24a6727a3a/sierralocal-0.1.0-py3-none-any.whl (8.0MB)
    100% |████████████████████████████████| 8.0MB 255kB/s 
   Installing collected packages: sierralocal
   Successfully installed sierralocal-0.1.0
   ```
   
   However, attempting to run the program will throw this exception:
   ```console
   art@Jesry:~$ sierralocal RT.fa
   searching path /usr/local/lib/python3.6/dist-packages/sierralocal/data/HIVDB*.xml
   Error: could not find local copy of HIVDB XML, attempting download...
   Unable to update HIVDB XML. Try manually downloading the HIVdb ASI2.
   [Errno 13] Permission denied: '/usr/local/lib/python3.6/dist-packages/sierralocal/data/HIVDB_8.7.9e470b87.xml'
   /usr/local/lib/python3.6/dist-packages/sierralocal/data/apobec.tsv
   Error: could not retrieve APOBEC DRM data
   Unable to update APOBEC DRMs. Try manually downloading the APOBEC DRM TSV into data/apobec.tsv
   Traceback (most recent call last):
     [...omitted...]
   FileNotFoundError: [Errno 2] No such file or directory: 'None'
   ```
   This occurs because Python does not have permission to write the HIVdb XML files to `/usr/local/lib`:
   ```console
   art@Jesry:~$ ls -l /usr/local/lib/python3.6/dist-packages/sierralocal
   total 100
   drwxr-sr-x 2 root staff  4096 Nov 15 10:59 bin
   drwxr-sr-x 2 root staff  4096 Nov 15 10:59 data
   -rw-r--r-- 1 root staff 13386 Nov 15 10:59 hivdb.py
   ```
   To rectify this, we have to elevate the read/write/execute privileges to everyone:
   ```console
   art@Jesry:~$ sudo chmod 777 /usr/local/lib/python3.6/dist-packages/sierralocal/data
   art@Jesry:~$ ls -l /usr/local/lib/python3.6/dist-packages/sierralocal
   total 100
   drwxr-sr-x 2 root staff  4096 Nov 15 10:59 bin
   drwxrwsrwx 2 root staff  4096 Nov 15 10:59 data
   -rw-r--r-- 1 root staff 13386 Nov 15 10:59 hivdb.py
   ```
   The *sierra-local* `setup.py` script does this automatically, which is why it is the recommended method if you have sudo privileges.
   Now the program should run properly:
   ```console
   art@Jesry:~$ sierralocal RT.fa
   searching path /usr/local/lib/python3.6/dist-packages/sierralocal/data/HIVDB*.xml
   Error: could not find local copy of HIVDB XML, attempting download...
   Updated HIVDB XML from https://hivdb.stanford.edu/assets/media/HIVDB_8.7.9e470b87.xml into /usr/local/lib/python3.6/dist-packages/sierralocal/data/HIVDB_8.7.9e470b87.xml
   /usr/local/lib/python3.6/dist-packages/sierralocal/data/apobec.tsv
   Error: could not retrieve APOBEC DRM data
   Updated APOBEC DRMs from https://hivdb.stanford.edu/assets/media/apobec-drms.5b7e1215.tsv into /usr/local/lib/python3.6/dist-packages/sierralocal/data/apobec.tsv
   HIVdb version 8.7
   Found NucAmino binary /usr/local/lib/python3.6/dist-packages/sierralocal/bin/nucamino-linux-amd64
   Aligned RT.fa
   100 sequences found in file RT.fa.
   Writing JSON to file git/sierra-local/RT_results.json
   Time elapsed: 8.1593 seconds (18.869 it/s)
   ```

## For macOS users


