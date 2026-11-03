[CmdletBinding()]
param(
    [ValidateRange(1, 6)]
    [int]$StartStep = 1,

    [ValidateRange(1, 6)]
    [int]$EndStep = 6,

    [switch]$DryRun,

    [switch]$CheckOnly,

    [switch]$RetryIncomplete,

    [switch]$ContinueOnError,

    [switch]$SkipPreflight,

    [string]$BashPath = "",

    [AllowEmptyString()]
    [string]$CondaEnv = "bioinfo_win"
)

$ErrorActionPreference = "Stop"
$script:ProjectRoot = $PSScriptRoot
$script:StateDirectory = Join-Path $script:ProjectRoot "_pipeline_state"
$script:LogDirectory = Join-Path $script:ProjectRoot "_pipeline_logs"
$script:MainLogPath = Join-Path $script:LogDirectory "pipeline_master.log"
$script:RunStartedAt = Get-Date

function Write-PipelineMessage {
    param(
        [string]$Message,
        [string]$Color = "Gray"
    )

    Write-Host $Message -ForegroundColor $Color
    if (-not $CheckOnly -and -not $DryRun -and $null -ne $script:MainLogPath) {
        $timestamp = (Get-Date).ToString("yyyy-MM-dd HH:mm:ss")
        Add-Content -LiteralPath $script:MainLogPath -Value "[$timestamp] $Message" -Encoding UTF8
    }
}

function Convert-ToBashSingleQuoted {
    param([string]$Value)

    return "'" + $Value.Replace("'", "'""'""'") + "'"
}

function Resolve-BashExecutable {
    param([string]$RequestedPath)

    $candidates = New-Object System.Collections.Generic.List[string]
    if (-not [string]::IsNullOrWhiteSpace($RequestedPath)) {
        $candidates.Add($RequestedPath)
    }

    $gitCommand = Get-Command git.exe -ErrorAction SilentlyContinue
    if ($null -ne $gitCommand) {
        $gitDirectory = Split-Path (Split-Path $gitCommand.Source -Parent) -Parent
        $candidates.Add((Join-Path $gitDirectory "bin\bash.exe"))
        $candidates.Add((Join-Path $gitDirectory "usr\bin\bash.exe"))
    }

    $candidates.Add("D:\Git\bin\bash.exe")
    $candidates.Add("C:\Program Files\Git\bin\bash.exe")
    $candidates.Add("C:\Program Files (x86)\Git\bin\bash.exe")

    foreach ($candidate in $candidates) {
        if (-not [string]::IsNullOrWhiteSpace($candidate) -and (Test-Path -LiteralPath $candidate -PathType Leaf)) {
            return (Resolve-Path -LiteralPath $candidate).Path
        }
    }

    throw "没有找到 Git Bash。请安装 Git for Windows，或使用 -BashPath 指定 bash.exe。"
}

function Get-PipelineSpecies {
    $listPath = Join-Path $script:ProjectRoot "folder_list.txt"
    if (Test-Path -LiteralPath $listPath -PathType Leaf) {
        return @(
            Get-Content -LiteralPath $listPath |
                ForEach-Object { $_.Trim() } |
                Where-Object { -not [string]::IsNullOrWhiteSpace($_) } |
                Select-Object -Unique
        )
    }

    return @(
        Get-ChildItem -LiteralPath $script:ProjectRoot -Directory |
            Where-Object {
                $_.Name -notlike "Auto_*" -and
                $_.Name -notlike "_pipeline_*" -and
                $_.Name -notlike ".*"
            } |
            ForEach-Object { $_.Name } |
            Sort-Object
    )
}

function Get-SpeciesMap {
    $mapPath = Join-Path $script:ProjectRoot "species_map.tsv"
    if (-not (Test-Path -LiteralPath $mapPath -PathType Leaf)) {
        return @()
    }

    $entries = New-Object System.Collections.Generic.List[object]
    foreach ($line in (Get-Content -LiteralPath $mapPath -Encoding UTF8)) {
        if ([string]::IsNullOrWhiteSpace($line)) {
            continue
        }
        $columns = $line -split "`t"
        if ($columns.Count -lt 6) {
            continue
        }
        $sourceDirectory = Join-Path $script:ProjectRoot $columns[1]
        $entries.Add([pscustomobject]@{
            SpeciesId = $columns[0]
            SourceDirectory = $sourceDirectory
            PepName = $columns[2]
            CdsName = $columns[3]
            GffName = $columns[4]
            LensName = $columns[5]
        })
    }

    return $entries.ToArray()
}

function Get-DiscoveredSpeciesInputs {
    $entries = New-Object System.Collections.Generic.List[object]
    $seen = New-Object System.Collections.Generic.HashSet[string]

    foreach ($directory in (Get-ChildItem -LiteralPath $script:ProjectRoot -Directory)) {
        if ($directory.Name -like "Auto_*" -or
            $directory.Name -like "_pipeline_*" -or
            $directory.Name -like "_script_backup*" -or
            $directory.Name -like "_wgdi*" -or
            $directory.Name -like "_invalid*" -or
            $directory.Name -like ".*") {
            continue
        }

        $pepFiles = @(Get-ChildItem -LiteralPath $directory.FullName -File -Filter "*.pep" -ErrorAction SilentlyContinue)
        if ($pepFiles.Count -ne 1) {
            continue
        }
        $speciesId = [System.IO.Path]::GetFileNameWithoutExtension($pepFiles[0].Name)
        $requiredFiles = @(
            $pepFiles[0].FullName,
            (Join-Path $directory.FullName ($speciesId + ".cds")),
            (Join-Path $directory.FullName ($speciesId + ".new.gff")),
            (Join-Path $directory.FullName ($speciesId + ".lens"))
        )
        if (@($requiredFiles | Where-Object { Test-Path -LiteralPath $_ -PathType Leaf }).Count -ne 4) {
            continue
        }
        if ($seen.Add($speciesId)) {
            $entries.Add([pscustomobject]@{
                SpeciesId = $speciesId
                SourceDirectory = $directory.FullName
                PepName = $pepFiles[0].Name
                CdsName = $speciesId + ".cds"
                GffName = $speciesId + ".new.gff"
                LensName = $speciesId + ".lens"
            })
        }
    }

    return $entries.ToArray()
}

function Get-ComparisonFolders {
    param([string[]]$Species)

    $folders = New-Object System.Collections.Generic.List[object]
    foreach ($baseName in $Species) {
        $baseDirectory = Join-Path $script:ProjectRoot ("Auto_wgdi\" + $baseName + "_base")
        foreach ($otherName in $Species) {
            $folderName = $baseName + "_" + $otherName
            $folders.Add([pscustomobject]@{
                BaseName = $baseName
                OtherName = $otherName
                IsSame = ($baseName -eq $otherName)
                Name = $folderName
                Path = Join-Path $baseDirectory $folderName
            })
        }
    }

    return $folders.ToArray()
}

function New-ValidationResult {
    param(
        [string[]]$Errors,
        [string]$Summary
    )

    $actualErrors = @($Errors | Where-Object { -not [string]::IsNullOrWhiteSpace($_) })
    return [pscustomobject]@{
        IsValid = ($actualErrors.Count -eq 0)
        Errors = $actualErrors
        Summary = $Summary
    }
}

function Test-Step01 {
    $errors = New-Object System.Collections.Generic.List[string]
    $speciesMap = @(Get-SpeciesMap)
    $discoveredSpecies = @(Get-DiscoveredSpeciesInputs)

    $mappedIds = @($speciesMap | ForEach-Object { $_.SpeciesId } | Sort-Object)
    $discoveredIds = @($discoveredSpecies | ForEach-Object { $_.SpeciesId } | Sort-Object)
    if (($mappedIds -join "`n") -ne ($discoveredIds -join "`n")) {
        $errors.Add(
            "完整物种清单已变化：当前输入为 [$($discoveredIds -join ', ')]，" +
            "species_map.tsv 为 [$($mappedIds -join ', ')]。需要重新执行第 1 步。"
        )
    }

    if ($speciesMap.Count -eq 0) {
        $errors.Add("species_map.tsv 不存在或没有完整物种记录。")
        $summary = if ($discoveredSpecies.Count -gt 0) {
            "待生成物种映射（当前完整输入 $($discoveredSpecies.Count) 个）"
        }
        else {
            "没有完整物种"
        }
        return (New-ValidationResult $errors.ToArray() $summary)
    }

    foreach ($entry in $speciesMap) {
        foreach ($fileName in @($entry.PepName, $entry.CdsName, $entry.GffName, $entry.LensName)) {
            $filePath = Join-Path $entry.SourceDirectory $fileName
            if (-not (Test-Path -LiteralPath $filePath -PathType Leaf)) {
                $errors.Add("缺少输入文件：$filePath")
            }
        }
    }

    $foundCount = 0
    foreach ($baseSpecies in $speciesMap) {
        foreach ($targetSpecies in $speciesMap) {
            $blastPath = Join-Path $script:ProjectRoot (
                "Auto_pep\" + $baseSpecies.SpeciesId + "_pep\" +
                $baseSpecies.SpeciesId + "_" + $targetSpecies.SpeciesId + ".blast"
            )
            $movedBlastPath = Join-Path $script:ProjectRoot (
                "Auto_wgdi\" + $baseSpecies.SpeciesId + "_base\" +
                $baseSpecies.SpeciesId + "_" + $targetSpecies.SpeciesId + "\" +
                $baseSpecies.SpeciesId + "_" + $targetSpecies.SpeciesId + ".blast"
            )
            if ((Test-Path -LiteralPath $blastPath -PathType Leaf) -or
                (Test-Path -LiteralPath $movedBlastPath -PathType Leaf)) {
                $foundCount++
            }
            else {
                $errors.Add("缺少 BLAST 结果：$blastPath（也未在 WGDI 目录找到）")
            }
        }
    }

    $expectedCount = $speciesMap.Count * $speciesMap.Count
    $summary = "BLAST 文件 $foundCount/$expectedCount"
    return (New-ValidationResult $errors.ToArray() $summary)
}

function Test-Step02 {
    $errors = New-Object System.Collections.Generic.List[string]
    $species = @(Get-PipelineSpecies)
    if ($species.Count -eq 0) {
        $errors.Add("没有找到物种目录。")
        return (New-ValidationResult $errors.ToArray() "没有物种")
    }

    $comparisonFolders = @(Get-ComparisonFolders $species)
    foreach ($comparison in $comparisonFolders) {
        if (-not (Test-Path -LiteralPath $comparison.Path -PathType Container)) {
            $errors.Add("缺少比较文件夹：$($comparison.Path)")
        }
    }

    $summary = "比较文件夹 $($comparisonFolders.Count) 个"
    return (New-ValidationResult $errors.ToArray() $summary)
}

function Test-Step03 {
    $errors = New-Object System.Collections.Generic.List[string]
    $species = @(Get-PipelineSpecies)
    if ($species.Count -eq 0) {
        $errors.Add("没有找到物种目录。")
        return (New-ValidationResult $errors.ToArray() "没有物种")
    }

    $comparisonFolders = @(Get-ComparisonFolders $species)
    $found = 0
    foreach ($comparison in $comparisonFolders) {
        $blastPath = Join-Path $comparison.Path ($comparison.Name + ".blast")
        if (Test-Path -LiteralPath $blastPath -PathType Leaf) {
            $found++
        }
        else {
            $errors.Add("缺少已归位的 BLAST 文件：$blastPath")
        }
    }

    $summary = "归位 BLAST $found/$($comparisonFolders.Count)"
    return (New-ValidationResult $errors.ToArray() $summary)
}

function Test-ConfigStep {
    param(
        [ValidateSet("Different", "Same")]
        [string]$Mode
    )

    $errors = New-Object System.Collections.Generic.List[string]
    $species = @(Get-PipelineSpecies)
    if ($species.Count -eq 0) {
        $errors.Add("没有找到物种目录。")
        return (New-ValidationResult $errors.ToArray() "没有物种")
    }

    $comparisonFolders = @(Get-ComparisonFolders $species)
    if ($Mode -eq "Different") {
        $comparisonFolders = @($comparisonFolders | Where-Object { -not $_.IsSame })
    }
    else {
        $comparisonFolders = @($comparisonFolders | Where-Object { $_.IsSame })
    }

    if ($comparisonFolders.Count -eq 0) {
        if ($Mode -eq "Different") {
            return (New-ValidationResult @() "不同物种 WGDI 配置 0/0（单物种自比对）")
        }
        $errors.Add("没有找到需要生成配置的同物种比较目录。")
        return (New-ValidationResult $errors.ToArray() "同物种 WGDI 配置 0/0")
    }

    $complete = 0
    foreach ($comparison in $comparisonFolders) {
        $runPath = Join-Path $comparison.Path "run.sh"
        $configPath = Join-Path $comparison.Path "total.conf"
        if ((Test-Path -LiteralPath $runPath -PathType Leaf) -and (Test-Path -LiteralPath $configPath -PathType Leaf)) {
            $complete++
        }
        else {
            $errors.Add("配置不完整：$($comparison.Path)")
        }
    }

    $summary = "WGDI 配置 $complete/$($comparisonFolders.Count)"
    return (New-ValidationResult $errors.ToArray() $summary)
}

function Test-Step06 {
    $errors = New-Object System.Collections.Generic.List[string]
    $species = @(Get-PipelineSpecies)
    if ($species.Count -eq 0) {
        $errors.Add("没有找到物种目录。")
        return (New-ValidationResult $errors.ToArray() "没有物种")
    }

    $comparisonFolders = @(Get-ComparisonFolders $species)
    $complete = 0
    foreach ($comparison in $comparisonFolders) {
        $prefix = if ($comparison.IsSame) {
            $comparison.BaseName
        }
        else {
            $comparison.Name
        }

        $requiredNames = @(
            "$prefix.png",
            "$prefix.icl",
            "$prefix.ks",
            "$prefix.bi.csv",
            "$prefix.bk.png",
            "$prefix.kp.png",
            "$prefix.pf.png",
            "run.sh",
            "total.conf"
        )
        if ($comparison.IsSame) {
            $requiredNames += "$prefix.bi.c.csv"
        }

        $missing = @(
            $requiredNames | Where-Object {
                -not (Test-Path -LiteralPath (Join-Path $comparison.Path $_) -PathType Leaf)
            }
        )
        if ($missing.Count -eq 0) {
            $complete++
        }
        else {
            $errors.Add("WGDI 结果不完整：$($comparison.Name)（缺少 $($missing -join ', ')）")
        }
    }

    $summary = "WGDI 结果 $complete/$($comparisonFolders.Count)"
    return (New-ValidationResult $errors.ToArray() $summary)
}

$script:Steps = @(
    [pscustomobject]@{
        Number = 1
        Script = "1run_blast.sh"
        Title = "建立 BLAST 数据库并执行蛋白比对"
        Validate = { Test-Step01 }
        Evidence = {
            $path = Join-Path $script:ProjectRoot "Auto_pep"
            if (Test-Path -LiteralPath $path -PathType Container) {
                return @(Get-ChildItem -LiteralPath $path -Recurse -File -ErrorAction SilentlyContinue |
                    Where-Object { $_.Name -like "*.db.*" -or $_.Name -like "*.blast" }).Count -gt 0
            }
            return $false
        }
    },
    [pscustomobject]@{
        Number = 2
        Script = "2run_pre_wgdi.sh"
        Title = "创建 WGDI 比较目录并复制输入"
        Validate = { Test-Step02 }
        Evidence = {
            $path = Join-Path $script:ProjectRoot "Auto_wgdi"
            if (Test-Path -LiteralPath $path -PathType Container) {
                return @(Get-ChildItem -LiteralPath $path -Directory -Filter "*_base" -ErrorAction SilentlyContinue).Count -gt 0
            }
            return $false
        }
    },
    [pscustomobject]@{
        Number = 3
        Script = "3blast_to_wgdi.sh"
        Title = "把 BLAST 结果移动到对应比较目录"
        Validate = { Test-Step03 }
        Evidence = {
            $path = Join-Path $script:ProjectRoot "Auto_wgdi"
            if (Test-Path -LiteralPath $path -PathType Container) {
                return @(Get-ChildItem -LiteralPath $path -Recurse -File -Filter "*.blast" -ErrorAction SilentlyContinue).Count -gt 0
            }
            return $false
        }
    },
    [pscustomobject]@{
        Number = 4
        Script = "4run_wgdi_different.sh"
        Title = "生成不同物种比较的 WGDI 配置"
        Validate = { Test-ConfigStep "Different" }
        Evidence = {
            $species = @(Get-PipelineSpecies)
            return @(Get-ComparisonFolders $species |
                Where-Object {
                    -not $_.IsSame -and
                    (Test-Path -LiteralPath (Join-Path $_.Path "run.sh") -PathType Leaf)
                }).Count -gt 0
        }
    },
    [pscustomobject]@{
        Number = 5
        Script = "5run_wgdi_same.sh"
        Title = "生成同物种比较的 WGDI 配置"
        Validate = { Test-ConfigStep "Same" }
        Evidence = {
            $species = @(Get-PipelineSpecies)
            return @(Get-ComparisonFolders $species |
                Where-Object {
                    $_.IsSame -and
                    (Test-Path -LiteralPath (Join-Path $_.Path "run.sh") -PathType Leaf)
                }).Count -gt 0
        }
    },
    [pscustomobject]@{
        Number = 6
        Script = "6run_wgdi.sh"
        Title = "执行全部 WGDI 重复程序"
        Validate = { Test-Step06 }
        Evidence = {
            $species = @(Get-PipelineSpecies)
            return @(Get-ComparisonFolders $species |
                Where-Object {
                    (Test-Path -LiteralPath (Join-Path $_.Path "total.conf") -PathType Leaf) -and
                    @(Get-ChildItem -LiteralPath $_.Path -File -Filter "*.icl" -ErrorAction SilentlyContinue).Count -gt 0
                }).Count -gt 0
        }
    }
)

function Get-StepStatePath {
    param([int]$StepNumber)
    return (Join-Path $script:StateDirectory ("step{0}.json" -f $StepNumber))
}

function Read-StepState {
    param([int]$StepNumber)

    $path = Get-StepStatePath $StepNumber
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        return $null
    }

    try {
        $json = [System.IO.File]::ReadAllText($path, [System.Text.Encoding]::UTF8)
        return ($json | ConvertFrom-Json)
    }
    catch {
        return [pscustomobject]@{
            Status = "corrupt"
            Error = $_.Exception.Message
        }
    }
}

function Save-StepState {
    param(
        [int]$StepNumber,
        [hashtable]$State
    )

    $path = Get-StepStatePath $StepNumber
    $temporaryPath = $path + "." + $PID + ".tmp"
    $json = $State | ConvertTo-Json -Depth 8
    $utf8WithoutBom = New-Object System.Text.UTF8Encoding($false)
    [System.IO.File]::WriteAllText($temporaryPath, $json, $utf8WithoutBom)
    Move-Item -LiteralPath $temporaryPath -Destination $path -Force
}

function Get-FileSha256 {
    param([string]$Path)

    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
        return ""
    }

    $stream = [System.IO.File]::OpenRead($Path)
    $sha256 = [System.Security.Cryptography.SHA256]::Create()
    try {
        $hashBytes = $sha256.ComputeHash($stream)
        return (($hashBytes | ForEach-Object { $_.ToString("x2") }) -join "").ToUpperInvariant()
    }
    finally {
        $sha256.Dispose()
        $stream.Dispose()
    }
}

function Get-StringSha256 {
    param([string]$Value)

    $sha256 = [System.Security.Cryptography.SHA256]::Create()
    try {
        $hashBytes = $sha256.ComputeHash([System.Text.Encoding]::UTF8.GetBytes($Value))
        return (($hashBytes | ForEach-Object { $_.ToString("x2") }) -join "").ToUpperInvariant()
    }
    finally {
        $sha256.Dispose()
    }
}

function Get-PipelineInputFingerprint {
    $parts = New-Object System.Collections.Generic.List[string]
    $parts.Add("pipeline-input-v1")
    $entries = @(Get-DiscoveredSpeciesInputs | Sort-Object SpeciesId)
    foreach ($entry in $entries) {
        $files = @(
            [pscustomobject]@{ Type = "PEP"; Name = $entry.PepName; Path = (Join-Path $entry.SourceDirectory $entry.PepName) },
            [pscustomobject]@{ Type = "CDS"; Name = $entry.CdsName; Path = (Join-Path $entry.SourceDirectory $entry.CdsName) },
            [pscustomobject]@{ Type = "GFF"; Name = $entry.GffName; Path = (Join-Path $entry.SourceDirectory $entry.GffName) },
            [pscustomobject]@{ Type = "LENS"; Name = $entry.LensName; Path = (Join-Path $entry.SourceDirectory $entry.LensName) }
        )
        foreach ($file in $files) {
            $fileHash = Get-FileSha256 $file.Path
            $parts.Add("$($entry.SpeciesId)|$($entry.SourceDirectory)|$($file.Type)|$($file.Name)|$fileHash")
        }
    }

    return Get-StringSha256 ($parts -join "`n")
}

function Get-ScriptHash {
    param([string]$ScriptName)
    $path = Join-Path $script:ProjectRoot $ScriptName
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        return ""
    }

    return Get-FileSha256 $path
}

function Test-ProcessAlive {
    param([int]$ProcessId)

    if ($ProcessId -le 0) {
        return $false
    }

    return ($null -ne (Get-Process -Id $ProcessId -ErrorAction SilentlyContinue))
}

function Get-BashPath {
    param([string]$WindowsPath)

    if ($WindowsPath -match "^([A-Za-z]):\\(.*)$") {
        $drive = $Matches[1].ToLowerInvariant()
        $rest = $Matches[2].Replace("\", "/")
        return "/$drive/$rest"
    }

    return $WindowsPath.Replace("\", "/")
}

function Get-CondaCommandPrefix {
    param(
        [string]$BashExecutable,
        [string]$EnvironmentName
    )

    if ([string]::IsNullOrWhiteSpace($EnvironmentName)) {
        return ""
    }

    $condaExecutable = (
        & $BashExecutable -lc "type -P conda" 2>$null |
            Select-Object -Last 1
    )
    if ([string]::IsNullOrWhiteSpace($condaExecutable)) {
        throw "没有在 Git Bash 中找到 conda。可以用 -CondaEnv '' 关闭自动激活。"
    }

    $condaExecutable = $condaExecutable.Trim()
    $scriptsDirectory = $condaExecutable.Substring(0, $condaExecutable.LastIndexOf("/"))
    $condaRoot = $scriptsDirectory.Substring(0, $scriptsDirectory.LastIndexOf("/"))
    $condaSh = $condaRoot + "/etc/profile.d/conda.sh"
    $quotedCondaSh = Convert-ToBashSingleQuoted $condaSh
    $quotedEnvironment = Convert-ToBashSingleQuoted $EnvironmentName

    return "source $quotedCondaSh && conda activate $quotedEnvironment && "
}

function Get-StepDecision {
    param(
        [object]$Step,
        [string]$CurrentHash,
        [string]$CurrentInputFingerprint
    )

    $state = Read-StepState $Step.Number
    if ($null -ne $state -and $state.Status -eq "completed") {
        if ([string]$state.InputFingerprint -ne $CurrentInputFingerprint) {
            return [pscustomobject]@{
                Action = "Run"
                Reason = "源数据内容指纹已变化，自动重建本步骤"
                State = $state
            }
        }

        if ($state.SourceHash -eq $CurrentHash) {
            $validation = & $Step.Validate
            if ($validation.IsValid) {
                return [pscustomobject]@{
                    Action = "SkipCompleted"
                    Reason = "已经完成，且脚本和产物检查均正常"
                    State = $state
                }
            }
            if (-not $RetryIncomplete) {
                return [pscustomobject]@{
                    Action = "SkipInvalidated"
                    Reason = "已有完成记录，但输入或产物已变化；默认不重跑"
                    State = $state
                }
            }
            return [pscustomobject]@{
                Action = "Run"
                Reason = "已完成记录的产物已失效，指定 -RetryIncomplete 后重跑"
                State = $state
            }
        }

        if (-not $RetryIncomplete) {
            return [pscustomobject]@{
                Action = "SkipChanged"
                Reason = "已经完成，但源脚本后来被修改；默认跳过"
                State = $state
            }
        }

        return [pscustomobject]@{
            Action = "Run"
            Reason = "源脚本已修改，且指定了 -RetryIncomplete"
            State = $state
        }
    }

    if ($null -ne $state) {
        if ($state.Status -eq "running" -and (Test-ProcessAlive ([int]$state.RunnerPid))) {
            return [pscustomobject]@{
                Action = "SkipActive"
                Reason = "检测到该步骤仍在运行"
                State = $state
            }
        }

        if (-not $RetryIncomplete) {
            return [pscustomobject]@{
                Action = "SkipStarted"
                Reason = "检测到启动记录，但尚未完成；默认跳过"
                State = $state
            }
        }

        return [pscustomobject]@{
            Action = "Run"
            Reason = "指定了 -RetryIncomplete，重新执行未完成步骤"
            State = $state
        }
    }

    $validation = & $Step.Validate
    if ($validation.IsValid) {
        return [pscustomobject]@{
            Action = "ReconcileCompleted"
            Reason = "未找到状态文件，但产物检查完整"
            State = $null
        }
    }

    $hasEvidence = [bool](& $Step.Evidence)
    if ($hasEvidence -and -not $RetryIncomplete) {
        return [pscustomobject]@{
            Action = "SkipStarted"
            Reason = "未找到状态文件，但检测到历史运行痕迹；默认跳过"
            State = $null
        }
    }

    return [pscustomobject]@{
        Action = "Run"
        Reason = if ($hasEvidence) { "检测到历史痕迹，指定 -RetryIncomplete 后重新执行" } else { "尚未启动" }
        State = $null
    }
}

function Test-PipelinePrerequisites {
    param([int]$StepNumber)

    $errors = New-Object System.Collections.Generic.List[string]

    if ($StepNumber -eq 1) {
        $species = @(
            Get-ChildItem -LiteralPath $script:ProjectRoot -Directory |
                Where-Object {
                    $_.Name -notlike "Auto_*" -and
                    $_.Name -notlike "_pipeline_*" -and
                    $_.Name -notlike ".*"
                }
        )

        if ($species.Count -eq 0) {
            $errors.Add("没有找到任何物种目录。")
        }

        $completeSpecies = 0
        $incompleteSpecies = 0
        foreach ($directory in $species) {
            $pepFiles = @(Get-ChildItem -LiteralPath $directory.FullName -File -Filter "*.pep" -ErrorAction SilentlyContinue)
            if ($pepFiles.Count -ne 1) {
                $incompleteSpecies++
                continue
            }

            $speciesId = [System.IO.Path]::GetFileNameWithoutExtension($pepFiles[0].Name)
            $requiredFiles = @(
                $pepFiles[0].FullName,
                (Join-Path $directory.FullName ($speciesId + ".cds")),
                (Join-Path $directory.FullName ($speciesId + ".new.gff")),
                (Join-Path $directory.FullName ($speciesId + ".lens"))
            )
            $allPresent = $true
            foreach ($requiredFile in $requiredFiles) {
                if (-not (Test-Path -LiteralPath $requiredFile -PathType Leaf)) {
                    $allPresent = $false
                    break
                }
            }
            if ($allPresent) {
                $completeSpecies++
            }
            else {
                $incompleteSpecies++
            }
        }

        Write-PipelineMessage "  物种文件检查：完整 $completeSpecies 个，不完整 $incompleteSpecies 个。" "DarkGray"
        if ($completeSpecies -lt 1) {
            $errors.Add("至少需要 1 个同时具备 .pep、.cds、.new.gff 和 .lens 的物种。")
        }

        foreach ($executable in @(
            "D:\ProgramFiles\bioinfo\bin\blastp.exe",
            "D:\ProgramFiles\bioinfo\bin\makeblastdb.exe"
        )) {
            if (-not (Test-Path -LiteralPath $executable -PathType Leaf)) {
                $errors.Add("缺少 BLAST 程序：$executable")
            }
        }
    }

    if ($StepNumber -eq 2 -and -not (Test-Path -LiteralPath (Join-Path $script:ProjectRoot "folder_list.txt") -PathType Leaf)) {
        $errors.Add("第 2 步需要第 1 步生成的 folder_list.txt。")
    }

    if ($StepNumber -eq 3 -and -not (Test-Path -LiteralPath (Join-Path $script:ProjectRoot "Auto_pep") -PathType Container)) {
        $errors.Add("第 3 步需要 Auto_pep 文件夹。")
    }

    if ($StepNumber -ge 4 -and $StepNumber -le 6) {
        if (-not (Test-Path -LiteralPath (Join-Path $script:ProjectRoot "folder_list.txt") -PathType Leaf)) {
            $errors.Add("第 $StepNumber 步需要 folder_list.txt。")
        }
    }

    if ($StepNumber -eq 6 -and -not $SkipPreflight) {
        try {
            $prefix = Get-CondaCommandPrefix $script:BashExecutable $CondaEnv
            $version = (
                & $script:BashExecutable -lc ($prefix + "wgdi --version") 2>&1 |
                    Out-String
            ).Trim()
            if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($version)) {
                $errors.Add("WGDI 不可用。请检查 -CondaEnv 或 conda 环境。")
            }
            else {
                Write-PipelineMessage "  WGDI 环境：$version" "DarkGray"
            }
        }
        catch {
            $errors.Add($_.Exception.Message)
        }
    }

    return $errors.ToArray()
}

function Invoke-PipelineStep {
    param(
        [object]$Step,
        [string]$CurrentHash,
        [string]$CurrentInputFingerprint
    )

    $previousState = Read-StepState $Step.Number
    $attempt = 1
    $startedAt = (Get-Date).ToString("o")
    if ($null -ne $previousState) {
        if ($null -ne $previousState.Attempt) {
            $attempt = [int]$previousState.Attempt + 1
        }
        if (-not [string]::IsNullOrWhiteSpace([string]$previousState.StartedAt)) {
            $startedAt = [string]$previousState.StartedAt
        }
    }

    Save-StepState $Step.Number @{
        Step = $Step.Number
        Script = $Step.Script
        Title = $Step.Title
        Status = "running"
        Started = $true
        StartedAt = $startedAt
        LastAttemptAt = (Get-Date).ToString("o")
        Attempt = $attempt
        RunnerPid = $PID
        SourceHash = $CurrentHash
        InputFingerprint = $CurrentInputFingerprint
        Log = $null
        ExitCode = $null
        CompletedAt = $null
        Validation = $null
        Errors = @()
    }

    $timestamp = (Get-Date).ToString("yyyyMMdd-HHmmss")
    $logPath = Join-Path $script:LogDirectory ("step{0}-{1}-{2}.log" -f $Step.Number, $Step.Script, $timestamp)
    $environmentPrefix = Get-CondaCommandPrefix $script:BashExecutable $CondaEnv
    $scriptArgument = "./" + $Step.Script
    $bashCommand = $environmentPrefix + "bash " + (Convert-ToBashSingleQuoted $scriptArgument)

    Write-PipelineMessage "  开始执行：$($Step.Script)" "Cyan"
    Write-PipelineMessage "  日志：$logPath" "DarkGray"
    Add-Content -LiteralPath $logPath -Value ("[{0}] START {1}" -f (Get-Date).ToString("o"), $Step.Script) -Encoding UTF8

    $exitCode = 1
    try {
        Push-Location $script:ProjectRoot
        try {
            & $script:BashExecutable -lc $bashCommand 2>&1 |
                Tee-Object -FilePath $logPath -Append
            $exitCode = $LASTEXITCODE
        }
        finally {
            Pop-Location
        }
    }
    catch {
        Add-Content -LiteralPath $logPath -Value ("[{0}] EXCEPTION {1}" -f (Get-Date).ToString("o"), $_.Exception.ToString()) -Encoding UTF8
        Write-PipelineMessage "  执行异常：$($_.Exception.Message)" "Red"
        $exitCode = 1
    }

    Add-Content -LiteralPath $logPath -Value ("[{0}] EXIT {1}" -f (Get-Date).ToString("o"), $exitCode) -Encoding UTF8
    $validation = & $Step.Validate
    $completed = ($exitCode -eq 0 -and $validation.IsValid)

    Save-StepState $Step.Number @{
        Step = $Step.Number
        Script = $Step.Script
        Title = $Step.Title
        Status = if ($completed) { "completed" } else { "failed" }
        Started = $true
        StartedAt = $startedAt
        LastAttemptAt = (Get-Date).ToString("o")
        Attempt = $attempt
        RunnerPid = $null
        SourceHash = $CurrentHash
        InputFingerprint = $CurrentInputFingerprint
        Log = $logPath
        ExitCode = $exitCode
        CompletedAt = if ($completed) { (Get-Date).ToString("o") } else { $null }
        Validation = $validation.Summary
        Errors = @($validation.Errors)
    }

    if ($completed) {
        Write-PipelineMessage "  完成：$($validation.Summary)" "Green"
        return $true
    }

    Write-PipelineMessage "  失败：退出码 $exitCode，$($validation.Summary)" "Red"
    foreach ($message in @($validation.Errors | Select-Object -First 12)) {
        Write-PipelineMessage "    - $message" "Red"
    }
    if ($validation.Errors.Count -gt 12) {
        Write-PipelineMessage "    - 其余 $($validation.Errors.Count - 12) 项请查看状态文件和日志。" "Red"
    }
    return $false
}

function Save-ReconciledState {
    param(
        [object]$Step,
        [string]$CurrentHash,
        [string]$CurrentInputFingerprint
    )

    $validation = & $Step.Validate
    Save-StepState $Step.Number @{
        Step = $Step.Number
        Script = $Step.Script
        Title = $Step.Title
        Status = "completed"
        Started = $false
        StartedAt = $null
        LastAttemptAt = (Get-Date).ToString("o")
        Attempt = 0
        RunnerPid = $null
        SourceHash = $CurrentHash
        InputFingerprint = $CurrentInputFingerprint
        Log = $null
        ExitCode = $null
        CompletedAt = (Get-Date).ToString("o")
        Validation = $validation.Summary
        Errors = @()
        ReconciledFromExistingOutput = $true
    }
}

function Get-AllStepStatus {
    foreach ($step in $script:Steps) {
        $validation = & $step.Validate
        $state = Read-StepState $step.Number
        [pscustomobject]@{
            Number = $step.Number
            Script = $step.Script
            Valid = $validation.IsValid
            Summary = $validation.Summary
            State = if ($null -eq $state) { "未记录" } else { $state.Status }
        }
    }
}

try {
    if ($StartStep -gt $EndStep) {
        throw "-StartStep 不能大于 -EndStep。"
    }

    $script:BashExecutable = Resolve-BashExecutable $BashPath
    foreach ($step in $script:Steps) {
        $scriptPath = Join-Path $script:ProjectRoot $step.Script
        if (-not (Test-Path -LiteralPath $scriptPath -PathType Leaf)) {
            throw "缺少分析脚本：$scriptPath"
        }
    }
    $currentInputFingerprint = Get-PipelineInputFingerprint

    if ($CheckOnly) {
        Write-PipelineMessage "只检查状态，不启动任何程序。" "Yellow"
        Write-PipelineMessage "项目目录：$script:ProjectRoot" "DarkGray"
        Write-PipelineMessage "Git Bash：$script:BashExecutable" "DarkGray"
        foreach ($step in $script:Steps) {
            $hash = Get-ScriptHash $step.Script
            $decision = Get-StepDecision $step $hash $currentInputFingerprint
            $validation = & $step.Validate
            $mark = if ($validation.IsValid) { "完整" } else { "不完整" }
            $color = if ($validation.IsValid) { "Green" } else { "Yellow" }
            Write-PipelineMessage ("[{0}/6] {1} | {2} | {3} | 状态判断：{4}" -f
                $step.Number, $step.Script, $mark, $validation.Summary, $decision.Reason) $color
        }
        exit 0
    }

    if ($DryRun) {
        Write-PipelineMessage "演练模式：只显示计划，不启动程序，也不写状态。" "Yellow"
        foreach ($step in $script:Steps | Where-Object { $_.Number -ge $StartStep -and $_.Number -le $EndStep }) {
            $hash = Get-ScriptHash $step.Script
            $decision = Get-StepDecision $step $hash $currentInputFingerprint
            Write-PipelineMessage ("[{0}/6] {1} -> {2}（{3}）" -f
                $step.Number, $step.Script, $decision.Action, $decision.Reason)
        }
        exit 0
    }

    New-Item -ItemType Directory -Path $script:StateDirectory -Force | Out-Null
    New-Item -ItemType Directory -Path $script:LogDirectory -Force | Out-Null
    Add-Content -LiteralPath $script:MainLogPath -Value ("`r`n===== RUN {0} | steps {1}-{2} =====" -f
        $script:RunStartedAt.ToString("o"), $StartStep, $EndStep) -Encoding UTF8

    $mutex = New-Object System.Threading.Mutex($false, "Local\HXC-WGDI-Pipeline")
    $lockTaken = $false
    try {
        try {
            $lockTaken = $mutex.WaitOne(0)
        }
        catch [System.Threading.AbandonedMutexException] {
            $lockTaken = $true
        }

        if (-not $lockTaken) {
            Write-PipelineMessage "已有另一个自动分析总控正在运行，本次不重复启动。" "Red"
            exit 2
        }

        Write-PipelineMessage "HXC 自动分析总控" "White"
        Write-PipelineMessage "步骤范围：$StartStep-$EndStep | 项目目录：$script:ProjectRoot" "DarkGray"
        if ($RetryIncomplete) {
            Write-PipelineMessage "已启用：重新执行未完成步骤。" "Yellow"
        }

        $runFailed = $false
        $blockedByDependency = $false

        for ($number = $StartStep; $number -le $EndStep; $number++) {
            $step = $script:Steps | Where-Object { $_.Number -eq $number }
            $hash = Get-ScriptHash $step.Script
            $decision = Get-StepDecision $step $hash $currentInputFingerprint
            Write-PipelineMessage ("`n[{0}/6] {1} - {2}" -f $step.Number, $step.Script, $step.Title) "White"
            Write-PipelineMessage "  状态判断：$($decision.Reason)" "DarkGray"

            switch ($decision.Action) {
                "SkipCompleted" {
                    $validation = & $step.Validate
                    Write-PipelineMessage "  已跳过：$($validation.Summary)" "Green"
                    continue
                }
                "SkipChanged" {
                    Write-PipelineMessage "  已跳过：源脚本在完成记录之后发生变化。若要重跑，请加 -RetryIncomplete。" "Yellow"
                    $blockedByDependency = $true
                    break
                }
                "SkipInvalidated" {
                    Write-PipelineMessage "  已跳过：已有完成记录，但当前输入或产物已变化。若要重建，请加 -RetryIncomplete。" "Yellow"
                    $blockedByDependency = $true
                    break
                }
                "SkipActive" {
                    Write-PipelineMessage "  已跳过：该步骤仍由另一个进程执行。" "Yellow"
                    $blockedByDependency = $true
                    break
                }
                "SkipStarted" {
                    Write-PipelineMessage "  已跳过：该步骤以前启动过，但未确认完成。" "Yellow"
                    Write-PipelineMessage "  如需续跑，请明确加 -RetryIncomplete；未经确认不会重复启动。" "Yellow"
                    $blockedByDependency = $true
                    break
                }
                "ReconcileCompleted" {
                    Save-ReconciledState $step $hash $currentInputFingerprint
                    $validation = & $step.Validate
                    Write-PipelineMessage "  已登记完成：$($validation.Summary)" "Green"
                    continue
                }
                "Run" {
                    if ($number -gt $StartStep) {
                        $priorUnfinished = @(
                            for ($prior = $StartStep; $prior -lt $number; $prior++) {
                                $priorStep = $script:Steps | Where-Object { $_.Number -eq $prior }
                                $priorValidation = & $priorStep.Validate
                                if (-not $priorValidation.IsValid) {
                                    $priorStep.Script
                                }
                            }
                        )
                        if ($priorUnfinished.Count -gt 0) {
                            Write-PipelineMessage "  阻止启动：前置步骤未确认完成（$($priorUnfinished -join ', ')）。" "Red"
                            $blockedByDependency = $true
                            break
                        }
                    }

                    if (-not $SkipPreflight) {
                        $prerequisiteErrors = @(Test-PipelinePrerequisites $number)
                        if ($prerequisiteErrors.Count -gt 0) {
                            Write-PipelineMessage "  前置条件未通过，未启动 $($step.Script)：" "Red"
                            foreach ($message in @($prerequisiteErrors | Select-Object -First 20)) {
                                Write-PipelineMessage "    - $message" "Red"
                            }
                            $blockedByDependency = $true
                            break
                        }
                    }

                    $stepSucceeded = Invoke-PipelineStep $step $hash $currentInputFingerprint
                    if (-not $stepSucceeded) {
                        $runFailed = $true
                        if (-not $ContinueOnError) {
                            break
                        }
                    }
                }
            }
        }

        Write-PipelineMessage "`n最终状态" "White"
        $allStatus = @(Get-AllStepStatus)
        foreach ($status in $allStatus) {
            $mark = if ($status.Valid) { "完成" } else { "未完成" }
            $color = if ($status.Valid) { "Green" } else { "Yellow" }
            Write-PipelineMessage ("  {0}. {1} | {2} | 记录：{3}" -f
                $status.Number, $status.Script, $mark, $status.State) $color
        }

        $allComplete = @($allStatus | Where-Object { -not $_.Valid }).Count -eq 0
        if ($allComplete) {
            Write-PipelineMessage "全部 6 个步骤及关键产物检查通过。" "Green"
            exit 0
        }

        if ($blockedByDependency) {
            Write-PipelineMessage "流程已在安全位置暂停，没有重复启动未完成程序。" "Yellow"
        }
        elseif ($runFailed) {
            Write-PipelineMessage "有步骤执行失败，请查看 _pipeline_logs 和 _pipeline_state。" "Red"
        }
        exit 4
    }
    finally {
        if ($lockTaken) {
            $mutex.ReleaseMutex()
        }
        $mutex.Dispose()
    }
}
catch {
    Write-PipelineMessage ("总控错误：{0}" -f $_.Exception.Message) "Red"
    if (-not [string]::IsNullOrWhiteSpace($_.InvocationInfo.PositionMessage)) {
        Write-PipelineMessage $_.InvocationInfo.PositionMessage "DarkRed"
    }
    exit 1
}
