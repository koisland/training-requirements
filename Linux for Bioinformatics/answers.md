# Linux for Bioinformatics

### Orientation to Linux Files and Navigation
* Q1: What is your home directory?
	* `/home/ubuntu`

* Q2: What is the output of the command?
	* `hello_world.txt`

* Q3: What is the output of each `ls` command?
	* `my_folder`
		* No output 
	* `my_folder2`
		* `hello_world.txt`

* Q4: What is the output of each?
	* `my_folder`
		* No output 
	* `my_folder2`
		* No output
	* `my_folder3`
		* `hello_world.txt`

### Scripting in the CLI
* Q5: What editor did you use and what was the command to save your file changes?
	* I used `nano` to edit `hello_world.txt`.
	* The command to save file changes is `ctrl + x`, then `y` to confirm changes, and finally `Enter` to write the filename.

### Set up a protected sudoer account and connect with it
* Q6: What is the error?
	* Permission denied (publickey). This means that the `sudouser` does not have an authorized public key for the ssh connection.

* Q7: What is the solution?
	* I found a solution from the article [How do I add new user accounts with SSH access to my Amazon EC2 Linux instance?](https://aws.amazon.com/premiumsupport/knowledge-center/new-user-accounts-linux-instance/) on aws.amazon.com. 
	* This required the following steps to solve Q6:
		1. Creating a `.ssh/authorized_keys` file in `/home/sudouser`
		2. Getting the public key from my private key with `ssh-keygen -y -f brn_linux_project.pem`.
		3. Adding the ssh-rsa public key to `.ssh/authorized_keys`.

### Docker
* Q8: What does the `sudo docker run` part of the command do? and what does the `salmon swim` part of the command do?
	* `sudo docker run` runs a specified docker image in a container.
	* According to `salmon -h`, `salmon swim` performs a super-secret operation. 

### Set up a non-sudo user account
* Q9: What is the output of this command?
	* `serveruser is not in the sudoers file. This incident will be reported.`

### Miniconda
* Q10: What is the output of `flask --version`?
	* `Python 3.8.12
	   Flask 2.0.2
	   Werkzeug 2.0.2`

### Mamba
* Q11: What is the output of `mamba -V`?
	* `conda 4.11.0`

### Conda environments
* Q12: What is the output of `which python`?
	* `/home/serveruser/miniconda3/envs/py27/bin/python`

* Q13: What is the output of `which python` now?
	* `/home/serveruser/miniconda3/bin/python`

### Installing `salmon`
* Q14: What is the output of `salmon -h`?

```
	salmon v1.7.0

	Usage:  salmon -h|--help or
	        salmon -v|--version or
	        salmon -c|--cite or
	        salmon [--no-version-check] <COMMAND> [-h | options]

	Commands:
	     index      : create a salmon index
	     quant      : quantify a sample
	     alevin     : single cell analysis
	     swim       : perform super-secret operation
	     quantmerge : merge multiple quantifications into a single file
``` 

### Simple RNA-Seq analysis with `salmon`

#### Part 1: Generating the transcriptome index
* Q15: What does the `-o athal.fa.gz` part of the command do?
	* The `-o` flag allows for an output filename.
* Q16: What is a `.gz` file?
	* A file with the `.gz` fil extenzion is a archive and compressed file using gzip.
	* https://en.wikipedia.org/wiki/Gzip
* Q17: What does the `zcat` command do?
	* `zcat` uncompresses files to standard output.
* Q18: What does the `head` command do?
	* `head` shows a specified number of lines from the beginning of a given data source.
* Q19: What does the number `100` signify in the command?
	* The number specifies the number of lines from the start of the file to print.
* Q20: What is `|` doing?
	* It 'pipes' or gives the output of one command to the next command to be used.
* Q21: What is a `.fa` file? What is the file format used for?
	* The file extension `.fa` is used for fasta files and serves as a de facto standard for storing sequence data.

#### Part 2: Quantify RNA-Seq data
* Q22: What format are the downloaded sequencing reads in?
	* The reads are in `.sra` format.
* Q23: What is the total size of the disk?
	* 7.7G is the total size of the disk listed as `/dev/xvda1`.
* Q24: How much space is remaining on the disk?
	* 1.2G is left on the disk.
* Q25: What went wrong?
	* The process ran out of storage space.
* Q26: What was your solution?
	* To compress the output `fastq` file using the `--gzip` argument for `fastq-dump`.
