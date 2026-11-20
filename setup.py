import setuptools

licenseDict = {"license_files": ["LICENSE"]}
if int(setuptools.__version__.split(".", 1)[0]) >= 77:
    try:
        import packaging.licenses
    except ImportError:
        licenseDict["license"] = "LGPLv3"
    else:
        licenseDict["license_expression"] = "LGPL-3.0-or-later"
else:
    licenseDict["license"] = "LGPLv3"
setuptools.setup(**licenseDict)
