# Copyright 2017 Grzegorz Firlik
#
# This file is part of PyDesc.
#
# PyDesc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# PyDesc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PyDesc.  If not, see <http://www.gnu.org/licenses/>.

class Param(object):

    """

    """

    def __init__(self, p_name, value, default_value=None):
        self.name = p_name
        self.value = value
        self.default_value = default_value

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str(self.name)


class Branch(object):

    def __init__(self, branch_name):
        self.name = branch_name

    def new_branch(self, branch_name):
        setattr(self, branch_name, Branch(self.name + '.' + branch_name))

    def del_branch(self, branch_name):
        delattr(self, branch_name)

    def set(self, param_name, value, default_value=None):

        fget = lambda self: self._get_p(param_name)
        fset = lambda self, value: self._set_p(
            param_name, value, default_value)
        fdel = lambda self: self._delete_p(param_name)

        setattr(self.__class__, param_name, property(fget, fset, fdel))
        setattr(self, '_' + param_name,
                Param(self.name + '.' + param_name, value, default_value))

    def set_default(self, param_name, value):
        self.set(param_name, value, value)

    def get(self, param_name):
        return self._get_p(param_name)

    def delete(self, param_name):
        current_default = getattr(self, '_' + param_name).default_value
        if current_default is None:
            self._delete_p(param_name)
        else:
            self.set(param_name, current_default, current_default)

    def _set_p(self, param_name, value, default_value):
        setattr(self, '_' + param_name,
                Param(self.name + '.' + param_name, value, default_value))

    def _get_p(self, param_name):
        return getattr(self, '_' + param_name).value

    def _delete_p(self, param_name):
        delattr(self, '_' + param_name)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return str(self.name)

    def _list_dict(self, slownik):
        for k, v in self.__dict__.items():
            if k is 'name':
                slownik.update({str(v): 'branch_sig'})
            elif k.startswith("_"):
                slownik.update({str(v): v.value})
            else:
                dummy = {}
                slownik.update(v._list_dict(dummy))
        return slownik


class ConfigManager(object):

    @staticmethod
    def new_branch(branch_name):
        setattr(ConfigManager, branch_name, Branch(
            ConfigManager.__name__ + '.' + branch_name))

    @staticmethod
    def del_branch(branch_name):
        delattr(ConfigManager, branch_name)

    @staticmethod
    def load_config(config_filename):
        try:
            with open(config_filename) as cofig_file_handler:
                exec(cofig_file_handler.read())
        except IOError as e:
            print("Configuration file read error: ", e)

    @staticmethod
    def _get_dict():
        slownik = {}
        for v in ConfigManager.__dict__.values():
            if str(v).startswith('ConfigManager'):
                v._list_dict(slownik)
        import operator
        return sorted(slownik.iteritems(), key=operator.itemgetter(0))

    @staticmethod
    def show_config(print_limit=0):
        sorted_slownik = ConfigManager._get_dict()
        print('+ConfigManager')
        for k, v in sorted_slownik:
            if v is 'branch_sig':
                print ('\t' * (k.count('.') - 1)), '+' + k.rsplit(".", 1)[1]
            else:
                try:
                    if len(v) > print_limit > 0:
                        try:
                            tmp = str(v[:print_limit])
                            print ('\t' * (k.count('.') - 1)), '-' + k.rsplit(".", 1)[
                                1] + ' = ' + tmp + ' ... <<' + str(len(v) - print_limit) + ' more>>'
                            continue
                        except TypeError:
                            from itertools import islice
                            print ('\t' * (k.count('.') - 1)), '-' + k.rsplit(".", 1)[1] + ' =', dict(
                                islice(v.iteritems(), print_limit)), '... <<' + str(len(v) - print_limit) + ' more>>'
                            continue
                except:
                    pass
                print ('\t' * (k.count('.') - 1)), '-' + \
                    k.rsplit(".", 1)[1] + ' = ' + str(v)

    @staticmethod
    def save_config(config_filename=None):
        if config_filename is None:
            import time
            config_filename = 'pydesc_config_' + time.strftime("%y%m%d-%H%M%S")
        sorted_slownik = ConfigManager._get_dict()
        try:
            with open(config_filename, 'w') as fh:
                for k, v in sorted_slownik:
                    if v is 'branch_sig':
                        tmp_a, tmp_b = k.rsplit(".", 1)
                        fh.write('%s.new_branch("%s")\n' % (tmp_a, tmp_b))
                    else:
                        tmp_a, tmp_b = str(k).rsplit(".", 1)
                        fh.write('%s.set("%s", %s)\n' % (tmp_a, tmp_b, str(v)))
        except IOError as e:
            print("Configuration file write error: ", e)
