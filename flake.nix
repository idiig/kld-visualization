{
  description = "Independent R plotting, Shiny and Shinylive Pages build environment";
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  inputs.flake-utils.url = "github:numtide/flake-utils";
  outputs = { nixpkgs, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        rEnv = pkgs.rWrapper.override {
          packages = with pkgs.rPackages; [ shiny ggrepel data_table igraph jsonlite httpuv shinylive ];
        };
        fonts = [ pkgs.noto-fonts-cjk-serif pkgs.noto-fonts-cjk-sans ];
      in {
        devShells.default = pkgs.mkShell {
          packages = [ rEnv pkgs.python3 pkgs.nodejs ] ++ fonts;
          FONTCONFIG_FILE = pkgs.makeFontsConf { fontDirectories = fonts; };
          KLD_PLOT_FONT = "Noto Serif CJK JP";
        };
      });
}
