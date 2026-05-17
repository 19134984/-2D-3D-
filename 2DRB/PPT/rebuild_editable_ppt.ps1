param(
    [string]$InputPptx = "",
    [string]$OutputPptx = "",
    [string]$LayoutJson = ".\PPT\_ocr_editable_layout.json"
)

$ErrorActionPreference = "Stop"
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8

function Resolve-FullPath([string]$PathText) {
    $executionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($PathText)
}

function Convert-HexToOleColor([string]$Hex) {
    $h = $Hex.TrimStart("#")
    $r = [Convert]::ToInt32($h.Substring(0, 2), 16)
    $g = [Convert]::ToInt32($h.Substring(2, 2), 16)
    $b = [Convert]::ToInt32($h.Substring(4, 2), 16)
    return [int]($r + 256 * $g + 65536 * $b)
}

function Add-EditableTextLayer($Slide, $Item, [bool]$AddMask, [bool]$TightMask) {
    $msoFalse = 0
    $msoTrue = -1
    $msoShapeRectangle = 1
    $msoTextOrientationHorizontal = 1

    if ($TightMask) {
        $padX = 1.0
        $padY = 0.5
    }
    else {
        $isStrongText = ([bool]$Item.bold) -or ([string]$Item.fg -eq "#0046A8") -or ([string]$Item.fg -eq "#FFFFFF")
        if ($isStrongText) {
            $padX = [Math]::Max(14.0, [double]$Item.height * 0.60)
            $padY = [Math]::Max(2.2, [double]$Item.height * 0.20)
        }
        else {
            $padX = [Math]::Max(7.0, [double]$Item.height * 0.35)
            $padY = [Math]::Max(1.8, [double]$Item.height * 0.16)
        }
    }
    $left = [double]$Item.left
    $top = [double]$Item.top
    $width = [Math]::Max(4.0, [double]$Item.width)
    $height = [Math]::Max(4.0, [double]$Item.height)

    if ($AddMask) {
        $mask = $Slide.Shapes.AddShape(
            $msoShapeRectangle,
            $left - $padX,
            $top - $padY,
            $width + 2 * $padX,
            $height + 2 * $padY
        )
        $mask.Name = "OCR mask"
        $mask.Fill.Visible = $msoTrue
        $mask.Fill.Solid()
        $mask.Fill.ForeColor.RGB = Convert-HexToOleColor ([string]$Item.bg)
        $mask.Line.Visible = $msoFalse
    }

    $box = $Slide.Shapes.AddTextbox(
        $msoTextOrientationHorizontal,
        $left - 0.2,
        $top - $height * 0.06,
        $width + [Math]::Max(4.0, $height * 0.25),
        $height * 1.18
    )
    $box.Name = "Editable OCR text"
    $box.Fill.Visible = $msoFalse
    $box.Line.Visible = $msoFalse
    $box.TextFrame2.MarginLeft = 0
    $box.TextFrame2.MarginRight = 0
    $box.TextFrame2.MarginTop = 0
    $box.TextFrame2.MarginBottom = 0
    $box.TextFrame2.WordWrap = $msoFalse
    try { $box.TextFrame2.VerticalAnchor = 3 } catch {}
    $range = $box.TextFrame2.TextRange
    $range.Text = [string]$Item.text
    $range.Font.Name = "Microsoft YaHei"
    $range.Font.NameFarEast = "Microsoft YaHei"
    $range.Font.Size = [double]$Item.fontSize
    $range.Font.Bold = $(if ([bool]$Item.bold) { $msoTrue } else { $msoFalse })
    $range.Font.Fill.ForeColor.RGB = Convert-HexToOleColor ([string]$Item.fg)
}

function Add-WhiteBase($Slide, [double]$Width, [double]$Height) {
    $msoFalse = 0
    $msoTrue = -1
    $msoShapeRectangle = 1
    $base = $Slide.Shapes.AddShape($msoShapeRectangle, 0, 0, $Width, $Height)
    $base.Name = "Editable white base"
    $base.Fill.Visible = $msoTrue
    $base.Fill.Solid()
    $base.Fill.ForeColor.RGB = Convert-HexToOleColor "#FFFFFF"
    $base.Line.Visible = $msoFalse
}

function Add-MovableAsset($Slide, $Asset, [string]$StageAssetDir) {
    $msoFalse = 0
    $msoTrue = -1
    $leaf = Split-Path -Leaf ([string]$Asset.file)
    $assetPath = Join-Path $StageAssetDir $leaf
    if (-not (Test-Path -LiteralPath $assetPath)) {
        return
    }
    $pic = $Slide.Shapes.AddPicture(
        $assetPath,
        $msoFalse,
        $msoTrue,
        [double]$Asset.left,
        [double]$Asset.top,
        [double]$Asset.width,
        [double]$Asset.height
    )
    $pic.Name = "Movable extracted image"
}

function Add-CleanupPatch($Slide, $Patch) {
    $msoFalse = 0
    $msoTrue = -1
    $msoShapeRectangle = 1
    $shape = $Slide.Shapes.AddShape(
        $msoShapeRectangle,
        [double]$Patch.left,
        [double]$Patch.top,
        [double]$Patch.width,
        [double]$Patch.height
    )
    $shape.Name = "Visual cleanup patch"
    $shape.Fill.Visible = $msoTrue
    $shape.Fill.Solid()
    $shape.Fill.ForeColor.RGB = Convert-HexToOleColor ([string]$Patch.bg)
    $shape.Line.Visible = $msoFalse
}

function Add-NativeShape($Slide, $ShapeDef) {
    $msoFalse = 0
    $msoTrue = -1
    $msoShapeRoundedRectangle = 5
    if ([string]$ShapeDef.type -eq "roundRect") {
        $shape = $Slide.Shapes.AddShape(
            $msoShapeRoundedRectangle,
            [double]$ShapeDef.left,
            [double]$ShapeDef.top,
            [double]$ShapeDef.width,
            [double]$ShapeDef.height
        )
        $shape.Name = "Editable native frame"
        $shape.Fill.Visible = $msoFalse
        $shape.Line.Visible = $msoTrue
        $shape.Line.ForeColor.RGB = Convert-HexToOleColor ([string]$ShapeDef.line)
        $shape.Line.Weight = [double]$ShapeDef.weight
    }
    elseif ([string]$ShapeDef.type -eq "line") {
        $shape = $Slide.Shapes.AddLine(
            [double]$ShapeDef.left,
            [double]$ShapeDef.top,
            [double]$ShapeDef.left + [double]$ShapeDef.width,
            [double]$ShapeDef.top + [double]$ShapeDef.height
        )
        $shape.Name = "Editable native line"
        $shape.Line.Visible = $msoTrue
        $shape.Line.ForeColor.RGB = Convert-HexToOleColor ([string]$ShapeDef.line)
        $shape.Line.Weight = [double]$ShapeDef.weight
    }
}

$pptDir = Resolve-FullPath ".\PPT"
if ([string]::IsNullOrWhiteSpace($InputPptx)) {
    $sourcePpt = Get-ChildItem -LiteralPath $pptDir -Filter "*.pptx" |
        Where-Object { $_.BaseName -notlike "*_*" } |
        Sort-Object LastWriteTime -Descending |
        Select-Object -First 1
    if ($null -eq $sourcePpt) {
        throw "No source pptx found in $pptDir"
    }
    $inputPath = $sourcePpt.FullName
}
else {
    $inputPath = Resolve-FullPath $InputPptx
}

if ([string]::IsNullOrWhiteSpace($OutputPptx)) {
    $suffixCodes = @(0x005F, 0x53EF, 0x7F16, 0x8F91, 0x91CD, 0x5EFA)
    $suffix = -join ($suffixCodes | ForEach-Object { [char]$_ })
    $outputPath = Join-Path (Split-Path -Parent $inputPath) (([IO.Path]::GetFileNameWithoutExtension($inputPath)) + $suffix + ".pptx")
}
else {
    $outputPath = Resolve-FullPath $OutputPptx
}
$layoutPath = Resolve-FullPath $LayoutJson
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$stageRoot = Join-Path $env:TEMP ("ppt-editable-rebuild-" + [Guid]::NewGuid().ToString("N"))
$stageInput = Join-Path $stageRoot "input.pptx"
$stageOutput = Join-Path $stageRoot "output.pptx"
$stagePreviewDir = Join-Path $stageRoot "preview"
$stageAssetDir = Join-Path $stageRoot "assets"

if (-not (Test-Path -LiteralPath $layoutPath)) {
    python (Join-Path $scriptDir "ocr_to_editable_layout.py")
}

$layout = Get-Content -LiteralPath $layoutPath -Encoding UTF8 -Raw | ConvertFrom-Json

New-Item -ItemType Directory -Path $stageRoot | Out-Null
New-Item -ItemType Directory -Path $stageAssetDir | Out-Null
Copy-Item -LiteralPath $inputPath -Destination $stageInput -Force
$assetSourceDir = Join-Path (Split-Path -Parent $layoutPath) "_editable_assets"
if (Test-Path -LiteralPath $assetSourceDir) {
    Get-ChildItem -LiteralPath $assetSourceDir -Filter "*.png" | ForEach-Object {
        Copy-Item -LiteralPath $_.FullName -Destination (Join-Path $stageAssetDir $_.Name) -Force
    }
}

$app = $null
$presentation = $null
try {
    $app = New-Object -ComObject PowerPoint.Application
    $app.Visible = -1
    $presentation = $app.Presentations.Open($stageInput, 0, 0, -1)
    $presentation.SaveAs($stageOutput, 24)

    foreach ($slideLayout in $layout.slides) {
        $page = [int]$slideLayout.page
        if ($page -gt 21) {
            continue
        }
        $slide = $presentation.Slides.Item($page)
        if ([string]$slideLayout.mode -eq "full_background") {
            $editableTextSet = @{}
            foreach ($text in $slideLayout.editableTexts) {
                $editableTextSet[[string]$text] = $true
            }
            foreach ($item in $slideLayout.items) {
                if ($editableTextSet.ContainsKey([string]$item.text)) {
                    Add-EditableTextLayer $slide $item $true $false
                }
            }
            Write-Host (
                "slide {0:00}: mode={1}, assets=0, text={2}" -f
                $page, $slideLayout.mode, $editableTextSet.Count
            )
            continue
        }
        if ([string]$slideLayout.mode -eq "white_asset") {
            Add-WhiteBase $slide ([double]$layout.slideWidth) ([double]$layout.slideHeight)
            foreach ($nativeShape in $slideLayout.nativeShapes) {
                Add-NativeShape $slide $nativeShape
            }
            foreach ($asset in $slideLayout.assets) {
                Add-MovableAsset $slide $asset $stageAssetDir
            }
            foreach ($patch in $slideLayout.cleanups) {
                Add-CleanupPatch $slide $patch
            }
        }
        $useTextMask = $true
        $tightTextMask = ([string]$slideLayout.mode -eq "white_asset")
        foreach ($item in $slideLayout.items) {
            Add-EditableTextLayer $slide $item $useTextMask $tightTextMask
        }
        Write-Host (
            "slide {0:00}: mode={1}, assets={2}, text={3}" -f
            $page, $slideLayout.mode, $slideLayout.assets.Count, $slideLayout.items.Count
        )
    }

    $presentation.Save()

    if (-not (Test-Path -LiteralPath $stagePreviewDir)) {
        New-Item -ItemType Directory -Path $stagePreviewDir | Out-Null
    }
    foreach ($page in @(1, 2, 10, 21, 22)) {
        $preview = Join-Path $stagePreviewDir ("slide{0:00}.png" -f $page)
        $presentation.Slides.Item($page).Export($preview, "PNG", 1672, 941) | Out-Null
    }

    Write-Host ("slide count: {0}" -f $presentation.Slides.Count)
}
finally {
    if ($presentation -ne $null) {
        $presentation.Close()
    }
    if ($app -ne $null) {
        $app.Quit()
    }
    [GC]::Collect()
    [GC]::WaitForPendingFinalizers()
}

Copy-Item -LiteralPath $stageOutput -Destination $outputPath -Force
$previewDir = Join-Path (Split-Path -Parent $outputPath) "_rebuild_preview"
if (-not (Test-Path -LiteralPath $previewDir)) {
    New-Item -ItemType Directory -Path $previewDir | Out-Null
}
Get-ChildItem -LiteralPath $stagePreviewDir -Filter "*.png" | ForEach-Object {
    Copy-Item -LiteralPath $_.FullName -Destination (Join-Path $previewDir $_.Name) -Force
}
Write-Host ("saved {0}" -f $outputPath)
