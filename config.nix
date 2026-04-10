{
  pkgs,
  config,
  ...
}:
{
  languages = {
    python = {
      enable = true;
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
        matplotlib
        scienceplots
        sphinx
      ];
    };
  };
}
