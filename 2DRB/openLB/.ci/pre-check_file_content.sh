#!/usr/bin/env nix-shell
#!nix-shell -p perl -i perl

use strict;
use warnings;

$ENV{LC_ALL} = "C";

sub check_whitespaces_and_tabs {
    my $src_endings = "c|cpp|cxx|h|hpp|hh";
    my @files = `git ls-files | grep -E "\\.($src_endings)\$"`;

    my $found_bad = 0;

    foreach my $file (@files) {
        chomp($file);
        my $content = `cat "$file"`;
        my @lines = split("\n", $content);

        for (my $i = 0; $i < scalar @lines; $i++) {
            my $line_number = $i + 1;
            my $line = $lines[$i];

            if ($line =~ /\t/) {
                print STDERR "Error: Indentation with Tab found in file $file at line $line_number.\n";
                $found_bad = 1;
            }

            if ($line =~ /[ \t]+$/) {
                print STDERR "Error: Trailing whitespace found in file $file at line $line_number.\n";
                $found_bad = 1;
            }
        }
    }

    if ($found_bad) {
        exit(1);
    }
}

sub check_nix_shebang_wrapping {
    my $found_bad = 0;

    open(my $ci, "<", ".gitlab-ci.yml") or die "Cannot read .gitlab-ci.yml";
    my $line_number = 0;
    while (my $line = <$ci>) {
        $line_number++;
        next unless $line =~ /nix develop \S+ --command (?:bash )?(\S+)/;
        my $script = $1;
        next unless -f $script;
        open(my $fh, "<", $script) or next;
        my $shebang = <$fh>;
        close($fh);
        if (defined $shebang && $shebang =~ m|^#!/usr/bin/env nix\s*$|) {
            print STDERR "Error: $script enters its own nix environment via shebang, ";
            print STDERR "but .gitlab-ci.yml line $line_number wraps it in nix develop.\n";
            $found_bad = 1;
        }
    }
    close($ci);

    if ($found_bad) {
        exit(1);
    }
}

sub check_no_compiled_binaries {
    my %magic = (
        "\x7fELF"          => "ELF",
        "\xcf\xfa\xed\xfe" => "Mach-O",
        "\xce\xfa\xed\xfe" => "Mach-O",
        "\xca\xfe\xba\xbe" => "Mach-O universal",
        "MZ"                => "PE",
    );

    my @files = `git ls-files`;

    my $found_bad = 0;

    foreach my $file (@files) {
        chomp($file);
        next unless -f $file;

        open(my $fh, "<:raw", $file) or next;
        my $head = "";
        read($fh, $head, 4);
        close($fh);

        foreach my $m (keys %magic) {
            next unless index($head, $m) == 0;
            print STDERR "Error: Compiled $magic{$m} binary committed: $file.\n";
            $found_bad = 1;
            last;
        }
    }

    if ($found_bad) {
        exit(1);
    }
}

# Do the work :-)

# Check whitespaces and tabs
check_whitespaces_and_tabs();

# Check that no nix shebang script is additionally wrapped in nix develop
check_nix_shebang_wrapping();

# Check that no compiled binary is committed
check_no_compiled_binaries();

# all OK
exit(0);
