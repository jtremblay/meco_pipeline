#!/usr/bin/env python

import re

class Mag:

    def __init__(self, name):
        if re.search("^\w[\w.-]*$", name):
            self._name = name
        else:
            raise Exception("Error: sample name \"" + name +
                "\" is invalid (should match [a-zA-Z0-9_][a-zA-Z0-9_.-]*)!")

        self._readsets = []

    def show(self):
        print("Mag -- name: " + self._name + ", readsets: " +
            ", ".join([readset.name for readset in self._readsets]))

    @property
    def name(self):
        return self._name

    @property
    def readsets(self):
        return self._readsets

    def readsets_by_name(self, name):
        return [readset for readset in self.readsets if readset.name == name]

    def add_magset(self, readset):
        if self.readsets_by_name(readset.name):
            raise Exception("Error: readset name \"" + readset.name +
                "\" already exists for sample \"" + self.name + "\"!")
        else:
            self.readsets.append(readset)
            readset._sample = self
