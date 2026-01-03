{
  pkgs,
  config,
  ...
}: {
  languages = {
    python = {
      enable = true;
      nixPackages = with pkgs."python${config.languages.python.version}Packages"; [hatchling];
    };
  };
}
