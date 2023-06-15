#!/usr/bin/env python

# Python Standard Modules
import configparser
import glob
import logging
import os
import re
import subprocess
import sys

log = logging.getLogger(__name__)

class Config(configparser.ConfigParser):

    def __init__(self):
        configparser.ConfigParser.__init__(self)
        # Do one single call to module avail when running pipeline wrapper.
        self._module_avail_output = subprocess.check_output(["bash", "-c", "module avail"], stderr=subprocess.STDOUT)
        #sys.stderr.write("module avail output: " + str(self._module_avail_output + "\n"))

    @property
    def filepath(self):
        return self._filepath

    def parse_files(self, config_files):
        self._filepath = config_files[-1].name

        # Make option names case sensitive
        self.optionxform = str
        for config_file in config_files:
            self.readfp(config_file)
        self.check_modules()

    # Check by a system call if all modules defined in config file are available
    def check_modules(self):
        modules = []

        # Retrieve all unique module version values in config file
        # assuming that all module key names start with "module_"
        for section in self.sections():
            for name, value in self.items(section):
                if re.search("^module_", name) and value not in modules:
                    modules.append(value)

        #sys.stderr.write("modules: " + str(modules) + "\n")
        log.info("Check modules...")
        for module in modules:
            # Bash shell must be invoked in order to find "module" cmd
            #module_show_output = subprocess.check_output(["bash", "-c", "module avail " + module], stderr=subprocess.STDOUT)
            # Actually, do not call module avail each time as some systems do not store modules in cache and can take 
            # very long to do multiple module av query each time we want to check if a module exists.
            if not (
                re.search(re.escape(module) + " ", self._module_avail_output.decode(), re.IGNORECASE) or 
                re.search(re.escape(module) + "\(default\)", self._module_avail_output.decode(), re.IGNORECASE) or 
                re.search(re.escape(module) + "\n", self._module_avail_output.decode(), re.IGNORECASE)
            ):
                sys.stderr.write("module : " + module + " was not found in module avail command. Please double check that desired module was spelled correctly.\n")
                raise Exception("Error in config file with " + module + "\n")
            else:
                log.info("Module " + module + " OK")
        log.info("Module check finished\n")

    # Retrieve param in config file with optional definition check and type validation
    # By default, parameter is required to be defined in the config file
    def param(self, section, option, required=True, type='string'):
        # Store original section for future error message, in case 'DEFAULT' section is used eventually
        original_section = section

        if not self.has_section(section):
            section = 'DEFAULT'

        if self.has_option(section, option):
            try:
                if type == 'int':
                    return self.getint(section, option)
                elif type == 'posint':
                    value = self.getint(section, option)
                    if value > 0:
                        return value
                    else:
                        raise Exception("Integer \"" + str(value) + "\" is not > 0!")
                elif type == 'float':
                    return self.getfloat(section, option)
                elif type == 'boolean':
                    return self.getboolean(section, option)
                elif type == 'filepath':
                    value = os.path.expandvars(self.get(section, option))
                    if os.path.isfile(value):
                        return value
                    else:
                        raise Exception("File path \"" + value + "\" does not exist or is not a valid regular file!")
                elif type == 'dirpath':
                    value = os.path.expandvars(self.get(section, option))
                    if os.path.isdir(value):
                        return value
                    else:
                        raise Exception("Directory path \"" + value + "\" does not exist or is not a valid directory!")
                elif type == 'prefixpath':
                    value = os.path.expandvars(self.get(section, option))
                    if glob.glob(value + "*"):
                        return value
                    else:
                        raise Exception("Prefix path \"" + value + "\" does not match any file!")
                elif type == 'list':
                    # Remove empty strings from list
                    return [x for x in self.get(section, option).split(",") if x]
                elif type == 'string':
                    return self.get(section, option)
                else:
                    raise Exception("Unknown parameter type '" + type + "'")
            except Exception as e:
                raise Exception("Error: parameter \"[" + section + "] " + option + "\" value \"" + self.get(section, option) + "\" is invalid!\n" + str(e))
        elif required:
            raise Exception("Error: parameter \"[" + original_section + "] " + option + "\" is not defined in config file!")
        else:
            return ""

# Global config object used throughout the whole pipeline
config = Config()
