# ----------------------------------------------------------------------------
# Copyright (c) 2025, AMeara.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

description = (
    "first draft"
)

setup(
    name="q2toy_1_cc",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    packages=find_packages(),
    author="AMeara",
    author_email="ameara@student.ethz.ch",
    description=description,
    url="https://example.com",
    entry_points={
        "qiime2.plugins": [
            "mytoy1="
            "mytoy1"
            ".plugin_setup:plugin"]
    },
    package_data={
        "mytoy1": ["citations.bib"],
        "mytoy1.tests": ["data/*"],
    },
    zip_safe=False,
)
