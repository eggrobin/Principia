$solutiondir = resolve-path $args[0]
$env:Path += ";$env:programfiles\Git\bin;$env:localappdata\GitHub\Portab~1\bin"
$newversion = (git describe --tags --always --dirty --abbrev=40 --long)
$headerpath = (join-path $solutiondir "base/version.hpp")
$cspath = (join-path $solutiondir "ksp_plugin_adapter/version.cs")
$date = (get-date).ToUniversalTime()

$generateversionheader = {
  $text = [string]::format(
      "#pragma once`n"                                `
          + "`n"                                      `
          + "namespace principia {{`n"                `
          + "namespace base {{`n"                     `
          + "`n"                                      `
          + "char const kBuildDate[] = `"{0:O}`";`n"  `
          + "char const kVersion[] =`n"               `
          + "    `"{1}`";`n"                          `
          + "`n"                                      `
          + "}}  // namespace base`n"                 `
          + "}}  // namespace principia`n",           `
      $date,
      $newversion)
  [system.io.file]::writealltext(
      $headerpath,
      $text,
      [system.text.encoding]::utf8)
  $cstext = [string]::format(
      "#pragma once`n"                                `
          + "`n"                                      `
          + "namespace principia {{`n"                `
          + "namespace ksp_plugin_adapter {{`n"       `
          + "`n"                                      `
          + "internal static class Version {{`n"      `
          + "  const string BuildDate = `"{0:O}`";`n" `
          + "  const string Version =`n"              `
          + "      `"{1}`";`n"                        `
          + "}}`n"                                    `
          + "`n"                                      `
          + "}}  // namespace base`n"                 `
          + "}}  // namespace principia`n",           `
      $date,
      $newversion)
  [system.io.file]::writealltext(
      $cspath,
      $cstext,
      [system.text.encoding]::utf8)
}

if (test-path -path $headerpath) {
  if ([system.io.file]::readalltext($headerpath) `
          -match '(?m)^\s+"([^"]+)";$.*') {
    $oldversion = $matches[1]
  }
  if ($oldversion.equals($newversion)) {
    echo "No change to git describe, leaving base/version.hpp untouched"
  } else {
    echo "Updating base/version.hpp, version is $newversion (was $oldversion)"
    &$generateversionheader
  }
} else {
  echo "Creating base/version.hpp, version is $newversion"
  &$generateversionheader
}
