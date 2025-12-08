#!/usr/bin/env python
# -*- encoding: utf-8 -*-


import select
import subprocess
import sys

__all__ = ["RunCommand", "cmd", "cmd_nowait"]


def _record_to_file(outfile, str):
    """
    Description: record str to the outfile
    Parameters:
        outfile: (str): the path of file
        str: (str): the string which to be recorded
    Returns:
        None
    """
    with open(outfile, "a") as fp:
        fp.write(str)


class RunCommand(object):
    """Provide some functions to run command or script
    This class will provide some method to run bash script or command
    by /bin/bash"""

    def __init__(self):
        self._cmd_str = None
        self._cwd = None
        self.communite_obj = None
        self._return_code = None
        self._stdout_list = []
        self._stderr_list = []
        self._start_run = False
        self.env = None

    @property
    def cmd_str(self):
        return self._cmd_str

    @cmd_str.setter
    def cmd_str(self, cmd):
        if isinstance(cmd, str):
            self._cmd_str = cmd
        else:
            raise TypeError("string/unicode is required")

    @property
    def cwd(self):
        return self._cwd

    @cwd.setter
    def cwd(self, cwd):
        if isinstance(cwd, str) or cwd is None:
            self._cwd = cwd
        else:
            raise TypeError("String is required")

    @property
    def return_code(self):
        return self._return_code

    @property
    def stdout(self):
        if self._stdout_list is None:
            return ""
        elif not isinstance(self._stdout_list, list):
            raise TypeError("STDOUT should be genrate from a list")
        else:
            return "".join(self._stdout_list)

    @property
    def stderr(self):
        if self._stderr_list is None:
            return ""
        elif not isinstance(self._stderr_list, list):
            raise TypeError("STDERR should be generate from a list")
        else:
            return "".join(self._stderr_list)

    @property
    def pid(self):
        if self.communite_obj is None:
            return None
        else:
            return self.communite_obj.pid

    @property
    def status(self):
        """The status of command or script
        Returns:
            0 (str), finished/running/notstart
        """
        status_finished = "finished"
        status_running = "running"
        status_not_start = "notstart"
        if self.start_run is False:
            return status_not_start
        elif self.communite_obj is None:
            return status_finished
        elif self.communite_obj.poll() is None:
            return status_running
        else:
            return status_finished

    def run(self):
        try:
            self.start_run = True
            self.communite_obj = subprocess.Popen(
                self.cmd_str,
                cwd=self.cwd,
                env=self.env,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                shell=True,
                executable="/bin/bash",
                # bufsize=1,  # line cache
            )
        except Exception as e:
            if self.communite_obj is not None:
                self.communite_obj.kill()
            self._return_code = -1
            self._stdout_list = []
            self._stderr_list = [str(e)]
            raise Exception(e)

    def terminate(self):
        """Terminate current command or script"""
        if self.communite_obj is not None:
            self.communite_obj.terminate()
            self._return_code = -1
            self.stderr = "%s\nTerminate by user" % self.stderr

    def process_stdout(self, line, stdout_list, display, outfile, check_error):
        stdout_list.append(line)
        self._stdout_list.append(line)
        if display:
            print(line, end="", file=sys.stdout)

        if outfile not in [None, ""]:
            _record_to_file(outfile, line)

        if check_error and line.strip() in check_error:
            self.terminate()

    def process_stderr(self, line, stderr_list, display, outfile, check_error):
        stderr_list.append(line)
        self._stderr_list.append(line)
        if display:
            print(line, end="", file=sys.stderr)

        if outfile not in [None, ""]:
            _record_to_file(outfile, line)

        if check_error and line.strip() in check_error:
            self.terminate()

    def get_std(
        self, tee=False, outfile="", display=True, check_error=None
    ):  # noqa: C901
        """Get stdout when running

        Returns:
            tuple: 1st is stdout string
                   2nd is stderr string
        """
        display = tee and display
        if self.start_run and self.communite_obj is not None:
            stdout_list = []
            stderr_list = []
            stdout = self.communite_obj.stdout
            stderr = self.communite_obj.stderr
            while True:
                try:
                    rlist, _, _ = select.select([stdout, stderr], [], [], 1)
                except ValueError:  # pipe closed
                    break

                if not rlist:
                    continue

                line_exists = False
                for stream in rlist:
                    if stream == stdout:
                        line = stdout.readline()
                        line = self.bytes_decode(line)
                        if line:
                            line_exists = True
                            self.process_stdout(
                                line, stdout_list, display, outfile, check_error
                            )
                    elif stream == stderr:
                        line = stderr.readline()
                        line = self.bytes_decode(line)
                        if line:
                            line_exists = True
                            self.process_stderr(
                                line, stderr_list, display, outfile, check_error
                            )

                if not line_exists:
                    break

            return "".join(stdout_list), "".join(stderr_list)

    def get_result(self):
        """Get the result of command/script
        Returns:
            tuple:
                1st is return code
                2nd is stdout
                3rd is stderr
            or
            None:
                command/script is not started/finished
        """
        # finished means the command/script is completed
        if self.status == "finished":
            if self._return_code is None:
                # self._return_code is not None means have been raised exception in funceion run()
                # if enter this choice means the command/script finished without exception
                # not means run successful or return 0
                self._return_code = self.communite_obj.poll()
                # put all stdout and stderr to self._stdout_list and self._stderr_list
                self.get_std()
            return (self.return_code, self.stdout, self.stderr)
        else:
            return None

    def bytes_decode(self, bytesstring):
        """change utf-8, gbk, big5, utf-16 to str"""
        encodinglist = [
            sys.stdout.encoding,
            sys.getfilesystemencoding(),
            "utf-8",
            "gbk",
            "big5",
            "utf-16",
        ]
        for encoding in encodinglist:
            if isinstance(encoding, bool):
                continue
            if encoding is None:
                continue
            try:
                content = bytesstring.decode(encoding)
                return content
            except UnicodeError:
                pass
        return bytesstring.decode(encoding=sys.stdout.encoding, errors="ignore")


def cmd(cmd, tee=True, outfile="", display=True, cwd=None, env=None, check_error=[]):
    """
    Description: Executes a command in a subprocess.
    Parameters:
        cmd: (str): bash-like command string.
        tee: (bool, optional): if True, simultaneously writes to stdout in real time while
             capturing output from the command. If not specified, defaults to
             True.
        outfile: (str, optional): file path, if not "" or None, it will record
                 all stdout and stderr into the file
        display: (bool, optional): if true, print command
        cwd (str, optional): if not None, change the working directory to cwd before executing
        check_error: (list, optional): stop running, if check_error in output/err
    Returns:
        return_code: (int): the return code of command line
        final_stdout: (str): the stdout of command line
        final_stderr: (str): the stderr of command line
    Example:
        >>> from cmd import cmd
        >>> cmd('ping -c 10.86.1.158')
        >>> cmd('ping -c 10.86.1.158', tee=False)
        >>> cmd('ping -c 10.86.1.158', tee=False, outfile='stdout.txt')
    """
    if display:
        print("[command] %s" % cmd)
    # Cannot display the output in Jenkins environment, must do sys.stdout.flush
    sys.stdout.flush()
    run_obj = RunCommand()
    run_obj.cmd_str = cmd
    run_obj.cwd = cwd
    run_obj.env = env
    run_obj.run()
    while True:
        run_obj.get_std(tee, outfile, display=display, check_error=check_error)
        # no display and record to log runtime
        if run_obj.status == "finished":
            break
    return_code, final_stdout, final_stderr = run_obj.get_result()
    return return_code, final_stdout, final_stderr


def cmd_nowait(cmd):
    """
    Description: Executes a command in a subprocess.
    Parameters:
        command: (str): bash-like command string.
    Returns:
        rub_obj: (instance of RunCommand): the instance of RunCommand
    Example:
        >>> from cmd import cmd_nowait
        >>> import os
        >>> import time
        >>> start_t = time.time()
        >>> if os.path.exists("1"):
        ...     os.remove("1")
        >>> cmd_nowait("sleep 5 && touch 1")
        >>> end_t = time.time()
        >>> assert(end_t-start_t < 3)
        >>> time.sleep(7)
        >>> assert(os.path.exists("1"))
        >>> os.remove("1")
    """
    print("[command no wait] %s" % cmd)
    # Cannot display the output in Jenkins environment, must do sys.stdout.flush
    sys.stdout.flush()
    run_obj = RunCommand()
    run_obj.cmd_str = cmd
    run_obj.run()
    return run_obj


if __name__ == "__main__":
    # cmd("echo hello world", tee=True)
    # result = cmd("sleep 5 && echo 1 && sleep 5 && echo 2", tee=False)
    # print("result:")
    # print(result)
    # result = cmd("sleep 5 && echo 1 && sleep 5 && echo 2", tee=True, outfile="test.log")
    # print("result:")
    # print(result)
    print(cmd("sleep 5 && echo 1"))
    cmd_nowait("sleep 5 && touch 1")
