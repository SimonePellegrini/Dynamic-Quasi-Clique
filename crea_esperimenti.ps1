# Verifica se sono stati passati nomi di dataset
if ($args.Count -eq 0) {
    Write-Host "Errore: Inserisci almeno un nome di dataset." -ForegroundColor Red
    Write-Host "Esempio: .\crea_esperimenti.ps1 dataset1 dataset2"
    exit
}

# Crea la cartella principale se non esiste
if (!(Test-Path "Esperimenti")) {
    New-Item -ItemType Directory -Path "Esperimenti" | Out-Null
}

foreach ($dataset in $args) {
    Write-Host "Elaborazione: $dataset" -ForegroundColor Cyan
    
    # Definisce i percorsi
    $pathStandard = "Esperimenti\$dataset\standard"
    $pathMixed = "Esperimenti\$dataset\mixed"
    
    # Crea le cartelle (l'opzione -Force crea anche le cartelle genitore)
    New-Item -ItemType Directory -Path $pathStandard -Force | Out-Null
    New-Item -ItemType Directory -Force -Path $pathMixed | Out-Null
    
    Write-Host "Cartelle create con successo per $dataset" -ForegroundColor Green
}