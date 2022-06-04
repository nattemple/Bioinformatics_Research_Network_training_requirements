# Answers to questions from "Linux for Bioinformatics"

Q1: 

In the command `ls`, what do `l` and `s` stand for? 

A: 

```
`ls stands for "**l**ist **s**torage"'
```

Q2: Change to the my_folder/ directory and then type ls. What is the output of this command?

A: 

```
hello_world.txt
```

Q3: What is the output of each ls command?

A:

The output of 'ls my_folder' is empty.
The output of 'ls my_folder2' is:

```
hello_world.txt
```

Q4. What is the output of each?

A:

The output of 'ls my_folder' is empty.
The output of 'ls my_folder2' is empty.
The output of 'ls my_folder3' is:

```
hello_world.txt
```

Q5. What editor did you use and what was the command to save your file changes?

A:

```
I used nano.  Within nano I pressed Cntl-S to save.
```

Q6. What is the error?

A:

```
SSH PUBLICKEY authentication failed.
```

Q7. What was the solution?

A:

`sudouser needed to have a .ssh/authorized_keys file.  I logged in as sudouser and created .ssh/authorized_keys with the following commands:`

```
$ mkdir .ssh
$ chmod 700 .ssh
$ touch .ssh/authorized_keys
$ chmod 600 .ssh/authorized_keys
```

I then logged out and returned to the ubuntu (root) user and copied authorized_keys to the sudouser account:

```
$ sudo cp .ssh/authorized_keys /home/sudouser/.ssh/
```

and changed its owner with:

```
$ sudo chown sudouser:sudouser /home/sudouser/.ssh/authorized_keys
```

Q8. what does the sudo docker run part of the command do? and what does the salmon swim part of the command do?

A:

```
The `sudo docker run` command creates a writeable container layer over the specified image (the `salmon swim` part).  The `salmon swim` part runs the command salmon swim (which could otherwise be run directly from the command line).
```

Q9. What is the output of this command?

A:

```
snap
```

Q10. What is the output of `flask --version`?

A:

```
Python 3.9.12
Flask 2.0.3
Werkzeug 2.0.3
```

Q11. What is the output of `mamba -V`

A: 

`conda 4.12.0`

Q12. What is the output of which python?

A:

`/home/serveruser/miniconda3/envs/py27/bin/python`

Q13. What is the output of which python now?

A:

`/home/serveruser/miniconda3/bin/python`

Q14. What is the output of `salmon -h

A:

```
salmon v1.4.0
```

Q15. What does the `-o athal.fa.gz` part of the command do?

A:

`Outputs to the file athal.fa.gz.`

Q16. What is a `.gz` file?

`A gzip file.`

Q17. What does the zcat command do?

A:

`The zcat uncompresses files to standard output.`

Q18. what does the head command do?

A:

`The head command displays the head (or first few lines) of a file.`

Q19. what does the number 100 signify in the command?

A:

`The top 100 lines.`

Q20. What is | doing? -- Hint using | in Linux is called "piping"

`It is piping the output of zcat athal.fa.gz to the input of head -n 100.`

Q21. What is a .fa file? What is this file format used for?

A:

`.fa files are FASTA files.  FASTA is a text-based format for representing nucleotide or peptide sequences.`

Q22. What format are the downloaded sequencing reads in?

A: 

`The downloaded files are in 'Sequence Read Archive' format.`

`The `prefetch` command produced these errors:`

```
$ prefetch SRR074122

Maximum file size download limit is 20,971,520KB

2022-05-29T03:20:38 prefetch.2.4.1 err: error unexpected while resolving tree within virtual file system module - failed to resolve accession 'SRR074122' - Obsolete software. See https://github.com/ncbi/sra-tools/wiki/Obsolete-software ( 406 )
2022-05-29T03:20:38 prefetch.2.4.1 err: error unexpected while resolving tree within virtual file system module - failed to resolve accession 'SRR074122' - Obsolete software. See https://github.com/ncbi/sra-tools/wiki/Obsolete-software ( 406 )
2022-05-29T03:20:38 prefetch.2.4.1 sys: error unexpected while resolving tree within virtual file system module - HTTP read failure: retrying...

Redirected!!!

2022-05-29T03:20:38 prefetch.2.4.1 err: name incorrect while evaluating path within network system module - Scheme is 'https'

Redirected!!!

2022-05-29T03:20:39 prefetch.2.4.1 err: name incorrect while evaluating path within network system module - Scheme is 'https'
2022-05-29T03:20:39 prefetch.2.4.1 sys: name incorrect while evaluating path within network system module - HTTP read failure: retrying...
2022-05-29T03:20:39 prefetch.2.4.1 err: path not found while resolving tree within virtual file system module - 'SRR074122' cannot be found.

```

`I instead downloaded the files using the link provided in the Readme.`

Q23. What is the total size of the disk?

A:

```
7.6 G
```

Q24. How much space is remaining on the disk?

A:

```
2.7 G
```

Q25. What went wrong?

A:

```
storage exhausted
```

Q26: What was your solution? 

A:

I modified the command to Write to gzip format:

```
$ fastq-dump --gzip ./SRR074122
```
