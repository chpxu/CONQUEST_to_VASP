{
  pkgs,
  config,
  ...
}:
{
  languages = {
    python = {
      enable = true;
      uv.enable = false;
      version = "313";
      nixPackages = with pkgs."python${config.languages.python.version}Packages"; [
        hatchling
        twine
        scipy
        numpy
        pytest
        flake8
        black
        pylint
        ase
        hatch
        matplotlib
        sphinx
        sphinx-rtd-theme
        keyring
      ];
    };
  };
}
