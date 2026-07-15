param(
    [string]$PandocPath = "$env:LOCALAPPDATA\Microsoft\WinGet\Packages\JohnMacFarlane.Pandoc_Microsoft.Winget.Source_8wekyb3d8bbwe\pandoc-3.10\pandoc.exe",
    [string]$ProjectRoot = '',
    [switch]$SkipImageRender
)

$ErrorActionPreference = 'Stop'
if ([string]::IsNullOrWhiteSpace($ProjectRoot)) {
    if ([string]::IsNullOrWhiteSpace($PSScriptRoot)) {
        throw 'ProjectRoot is required when the script is loaded from a UTF-8 ScriptBlock.'
    }
    $scriptDir = $PSScriptRoot
    $root = (Resolve-Path (Join-Path $scriptDir '..')).Path
}
else {
    $root = (Resolve-Path $ProjectRoot).Path
    $scriptDir = Join-Path $root '06_公式与质量验证'
}
$draftDir = Join-Path $root '03_推导草稿'
$figureDir = Join-Path $root '04_图表'
$finalDir = Join-Path $root '05_最终报告'
$tmpDir = Join-Path $root 'tmp\pdfs'
$imageDir = Join-Path $tmpDir 'images'
$finalMarkdown = Join-Path $finalDir 'TRT_magic参数与温度四阶误差_零基础教程.md'
$pdfInput = Join-Path $tmpDir 'TRT_report_pdf_input.md'
$finalPdf = Join-Path $finalDir 'TRT_magic参数与温度四阶误差_零基础教程.pdf'
$headerTex = Join-Path $scriptDir 'report_header.tex'
$utf8 = New-Object System.Text.UTF8Encoding($false)

New-Item -ItemType Directory -Force -Path $finalDir,$tmpDir,$imageDir | Out-Null

$chapters = @(
    '01_共同基础.md',
    '02_BGK到TRT.md',
    '03_边界magic.md',
    '04_温度四阶误差.md',
    '05_文献参数映射.md',
    '06_综合与练习.md'
)

$frontMatter = @'
---
title: "TRT magic 参数与温度四阶误差"
subtitle: "全新零基础版：从 population、碰撞与迁移到边界 3/16、D2Q5 的 1/6 和完整四阶消除"
date: "2026-07-14"
lang: zh-CN
---

# 阅读说明 {.unnumbered #reading-guide}

这是一份面向 **完全没有 LBM 基础** 的教程。报告从标量场、导数、格子、population、矩、碰撞和迁移开始，再逐步建立 BGK、TRT、Hénon 参数、边界表观滑移以及 D2Q5 体相四阶误差的完整关系。

核心结论只依据下列四篇本地 PDF：

1. Wang et al., *Lattice Boltzmann simulations of thermal convective flows in two dimensions*;
2. Contrino et al., *LB simulations of the thermally driven 2D square cavity at high Rayleigh numbers*;
3. Dubois and Lallemand, *Towards higher order lattice Boltzmann schemes*;
4. Ginzburg and d'Humières, *Multireflection boundary conditions for lattice Boltzmann models*.

本地第 2 篇 PDF 的文件名含 `[Luo2014]`，但论文第一作者是 **Dario Contrino**；本报告按论文实际作者书写。

## 先记住两条绝对不能混用的链 {.unnumbered}

- 受限的 D2Q9 流场 halfway bounce-back 边界链：$\Lambda^2=1/4\Longleftrightarrow H_{\mathrm{even}}H_{\mathrm{odd}}=3/16$。
- D2Q5 纯热体相四阶链：$H_\kappa H_\nu=1/6$ 只代表四阶各向同性；要同时消去两个指定四阶系数，还必须满足 Dubois Eq. (55) 中决定 $H_e$ 的关系。

这两条链的格子、场变量、误差对象和适用条件都不同。本报告不把局部体相系数消失宣称为整个热对流程序“全局四阶”，也不把受限平壁结论外推到任意边界。
'@

$parts = New-Object System.Collections.Generic.List[string]
$parts.Add($frontMatter.TrimEnd())

for ($i = 0; $i -lt $chapters.Count; $i++) {
    $path = Join-Path $draftDir $chapters[$i]
    $text = [System.IO.File]::ReadAllText($path, $utf8)
    $chapterId = $i + 1
    $text = [regex]::Replace(
        $text,
        '(?m)^# (.+)$',
        { param($m) '# ' + $m.Groups[1].Value + " {#chapter-$chapterId}" },
        1
    )
    $parts.Add($text.Trim())
}

$report = ($parts -join "`r`n`r`n") + "`r`n"
$report = $report.Replace('](01_共同基础.md)', '](#chapter-1)')
$report = $report.Replace('](02_BGK到TRT.md)', '](#chapter-2)')
$report = $report.Replace('](03_边界magic.md)', '](#chapter-3)')
$report = $report.Replace('](04_温度四阶误差.md)', '](#chapter-4)')
$report = $report.Replace('](05_文献参数映射.md)', '](#chapter-5)')
$report = $report.Replace('](06_综合与练习.md)', '](#chapter-6)')
[System.IO.File]::WriteAllText($finalMarkdown, $report, $utf8)

$pdfMarkdown = [regex]::Replace(
    $report,
    '\(\.\./04_图表/([^\)]+)\.svg\)',
    '(images/$1.png)'
)
[System.IO.File]::WriteAllText($pdfInput, $pdfMarkdown, $utf8)

if (-not $SkipImageRender) {
    $edge = 'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe'
    if (-not (Test-Path -LiteralPath $edge)) {
        throw "Microsoft Edge not found: $edge"
    }
    $edgeProfile = Join-Path $tmpDir 'edge-profile'
    New-Item -ItemType Directory -Force -Path $edgeProfile | Out-Null

    Get-ChildItem -LiteralPath $figureDir -Filter '*.svg' | ForEach-Object {
        $svgText = [System.IO.File]::ReadAllText($_.FullName, $utf8)
        $width = 1200
        $height = 800
        if ($svgText -match '<svg[^>]*\bwidth="([0-9.]+)"[^>]*\bheight="([0-9.]+)"') {
            $width = [math]::Ceiling([double]$Matches[1])
            $height = [math]::Ceiling([double]$Matches[2])
        }
        elseif ($svgText -match '<svg[^>]*\bviewBox="[^ ]+ [^ ]+ ([0-9.]+) ([0-9.]+)"') {
            $width = [math]::Ceiling([double]$Matches[1])
            $height = [math]::Ceiling([double]$Matches[2])
        }
        $png = Join-Path $imageDir ($_.BaseName + '.png')
        $uri = ([System.Uri]$_.FullName).AbsoluteUri
        $arguments = @(
            '--headless=new',
            '--disable-gpu',
            '--hide-scrollbars',
            '--no-first-run',
            '--force-device-scale-factor=2',
            "--user-data-dir=$edgeProfile",
            "--window-size=$width,$height",
            "--screenshot=$png",
            $uri
        )
        Start-Process -FilePath $edge -ArgumentList $arguments -WindowStyle Hidden -Wait | Out-Null
        for ($attempt = 0; $attempt -lt 50 -and -not (Test-Path -LiteralPath $png); $attempt++) {
            Start-Sleep -Milliseconds 100
        }
        if (-not (Test-Path -LiteralPath $png)) {
            throw "Failed to render SVG: $($_.FullName)"
        }
    }
}

if (-not (Test-Path -LiteralPath $PandocPath)) {
    throw "Pandoc not found: $PandocPath"
}

$texBin = 'D:\App\latex\texlive\2025\bin\windows'
$env:PATH = $texBin + ';' + $env:PATH
$resourcePath = $tmpDir + ';' + $root
$pandocArgs = @(
    $pdfInput,
    '--from=markdown+tex_math_dollars+tex_math_single_backslash+raw_tex+pipe_tables+task_lists',
    '--standalone',
    '--toc',
    '--toc-depth=3',
    '--top-level-division=chapter',
    '--pdf-engine=xelatex',
    "--resource-path=$resourcePath",
    "--include-in-header=$headerTex",
    '--variable=documentclass:ctexbook',
    '--variable=classoption:openany',
    '--variable=papersize:a4',
    '--variable=geometry:top=2.2cm,bottom=2.1cm,left=2.0cm,right=2.0cm',
    '--variable=CJKmainfont:Microsoft YaHei',
    '--variable=CJKsansfont:Microsoft YaHei',
    '--variable=CJKmonofont:SimSun',
    '--variable=mainfont:Times New Roman',
    '--variable=monofont:Consolas',
    '--variable=fontsize:11pt',
    '--variable=linestretch:1.18',
    '--variable=colorlinks:true',
    '--variable=linkcolor:blue',
    '--variable=urlcolor:blue',
    '--syntax-highlighting=tango',
    '--output',
    $finalPdf
)

& $PandocPath @pandocArgs
if ($LASTEXITCODE -ne 0) {
    throw "Pandoc/XeLaTeX failed with exit code $LASTEXITCODE"
}

Write-Output $finalMarkdown
Write-Output $finalPdf
