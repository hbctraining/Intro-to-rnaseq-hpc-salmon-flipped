# Workshop Schedule

> **NOTE:** The *Basic Data Skills* [Introduction to the command-line interface](https://hbctraining.github.io/Intro-to-Shell-flipped/schedule/) workshop is a prerequisite.
## Day 1

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 10:00 - 10:30 | [Workshop introduction](../lectures/Intro_to_workshop.pdf) | Radhika |
| 10:30 - 11:45 | [Introduction to Shell](../lessons/01_the_filesystem.md) | Radhika|
| 11:45 - 12:00 | Overview of self-learning materials and homework submission | Jihe |


### Self Learning #1

Before you start with the self-learning portion of the workshop, please check that **you are logged into O2** and **are working on a compute node** (i.e. your command prompt should have the word `compute` in it).

> If you are not logged into O2 or are not on a compute node, please follow the steps below as appropriate before you start with the self-learning lessons:
> 1. Log in using `ssh rc_trainingXX@o2.hms.harvard.edu` and enter your password (HMSXXcluster) (replace the "XX" in both the username and the password with the number you were assigned in class). 
> 2. Once you are on the login node, use `srun --pty -p interactive -t 0-2:00 --mem 1G /bin/bash` to get on a compute node.
> 3. Proceed with the self learning once your command prompt has the word `compute` in it.
> 4. If you log out between lessons (using the `exit` command twice), please follow points 1. and 2. above to log back in and get on a compute node when you restart with the self learning.

* [Wildcards and shortcuts in Shell](../lessons/02_wildcards_shortcuts.md)
* [Examining and creating files](../lessons/03_working_with_files.md)
* [Searching and redirection](../lessons/04_searching_files.md)
* [Shell scripts and variables in Shell](../lessons/05_shell-scripts_variable.md)

### Assignment #1
* All exercise questions from the self-learning lessons have been put together in a [text file](https://raw.githubusercontent.com/hbctraining/Intro-to-shell-flipped/master/homework/Day1_assignment.txt) (download for local access).
  * The text file can be opened with any text editor application (i.e. Notepad++, TextWrangler) on your local computer
* Add your solutions to the exercises in the downloaded .txt file and **upload the saved text file** to [Dropbox](https://www.dropbox.com/request/WG2XCQEmvMbw84H8s0PM) **day before the next class**.
* [Email us](mailto:hbctraining@hsph.harvard.edu) about questions regarding the homework that you need answered before the next class.
* Post questions that you would like to have reviewed in class [here](https://PollEv.com/hbctraining945).
* [Answer key](../homework/Day1_answer_key.txt)

***

## Day 2

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:45 | Self-learning lessons review | All |
| 10:45 - 12:00 | [Loops and automation](../lessons/06_loops_and_automation.md) | Meeta |

### Self Learning #2
* [Permissions and Environment Variables](../lessons/07_permissions_and_environment_variables.md)
* [Introduction to High-performance computing](../lessons/08_HPC_intro_and_terms.md)

### Assignment #2
* All exercise questions from the self-learning lessons have been put together in a [text file](https://raw.githubusercontent.com/hbctraining/Intro-to-shell-flipped/master/homework/Day2_assignment.txt) (download for local access).
  * The text file can be opened with any text editor application (i.e. Notepad++, TextWrangler) on your local computer
* Add your solutions to the exercises in the downloaded .txt file and **upload the saved text file** to [Dropbox](https://www.dropbox.com/request/6zPM8itMzGcTOgSxKQrc) **day before the next class**.
* [Email us](mailto:hbctraining@hsph.harvard.edu) about questions regarding the homework that you need answered before the next class.
* Post questions that you would like to have reviewed in class [here](https://PollEv.com/hbctraining945).
* [Answer key](../homework/Day2_answer_key.txt)

***

## Day 3

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:00 | Self-learning lessons review | All |
| 10:00 - 11:00 | [Introduction to the O2 cluster](../lectures/HPC_intro_O2_October2020.pdf)| Radhika |
| 11:00 - 11:30 | [Exercise](../activities/sbatch_exercise.md)|  |
| 11:30 - 11:45 | [Introduction to the O2 cluster - data storage](../lectures/HPC_intro_O2_October2020.pdf)| Radhika |
| 11:45 - 12:00 | [Wrap up](../lectures/shell-wrapup-nanocourse.pdf) | Radhika |

***

## Resources

Bash/Shell Cheat sheets:
* [http://fosswire.com/post/2007/08/unixlinux-command-cheat-sheet/](http://fosswire.com/post/2007/08/unixlinux-command-cheat-sheet/)
* [https://github.com/swcarpentry/boot-camps/blob/master/shell/shell_cheatsheet.md](https://github.com/swcarpentry/boot-camps/blob/master/shell/shell_cheatsheet.md)
* [tldr_ : Simplified version of the `man` pages (online and searchable)](https://tldr.ostera.io/)

## Advanced bash commands
If you are interested in learning some more advanced tools for working on the command-line, we encourage you to walk-through the materials linked below:

* [More fun with bash (advanced)](../lessons/extra_bash_tools.md)



***
*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *Some materials used in these lessons were derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
