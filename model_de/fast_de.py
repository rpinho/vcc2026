"""Optional local rank optimization is excluded because its upstream license is unspecified."""
def install_unchecked():
    raise RuntimeError("The unlicensed rank optimization is excluded; run without --fast")
def install():
    return install_unchecked()
