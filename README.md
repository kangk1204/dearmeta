# DearMeta

DearMeta는 Illumina EPIC GEO 메틸레이션 데이터셋을 다운로드하고 전처리/분석할 수 있는 커맨드라인 툴킷입니다. `dearmeta` 명령 한 번으로 GEO 메타데이터 수집부터 IDAT 다운로드, R 기반 분석 파이프라인 실행, 결과 보고서 생성까지 자동화합니다.

## 주요 기능
- GEO 시리즈 메타데이터 및 플랫폼 어노테이션 자동 수집
- IDAT 페어 다운로드 및 캐시 관리
- 전처리에 필요한 워크스페이스 디렉터리 구조 자동 생성
- R 스크립트(`scripts/analysis.R`)를 이용한 차등 메틸레이션 분석과 시각화
- 실행 설정(`runtime/configure.tsv`)과 로그 자동 보존

## 시스템 요구 사항
- Python 3.10 이상
- R 4.0 이상 권장 (analysis 스크립트 실행용)
- GEO API와 Supplementary Files 다운로드를 위한 인터넷 연결

R 패키지는 `scripts/install.R`로 설치할 수 있습니다.

```bash
Rscript scripts/install.R
```

## 설치

```bash
git clone https://github.com/kangk1204/dearmeta.git
cd dearmeta
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e .[dev]
```

## 빠른 시작
1. 분석하려는 GEO 시리즈 ID(예: `GSE123456`)를 준비합니다.
2. 필수 R 패키지를 설치합니다.
3. CLI를 실행합니다.

```bash
dearmeta download --gse GSE123456
# runtime/configure.tsv 파일의 dear_group 컬럼을 채운 뒤
dearmeta analysis --gse GSE123456
```

명령이 완료되면 각 시리즈별 디렉터리 내에 다음과 같은 구조가 만들어집니다.

```
GSE123456/
├── 01_download/        # 원본 IDAT 및 메타데이터
├── 02_preprocess/      # 전처리 산출물
├── 03_analysis/        # 분석 결과
├── 04_figures/         # 정적 플롯
├── 05_interactive/     # HTML 보고서
└── runtime/            # configure.tsv, 로그, 실행 설정
```

> **참고:** 실제 GEO 데이터와 분석 산출물은 용량이 매우 크므로 Git 저장소에는 포함하지 않습니다. `.gitignore`에 의해 `GSE*/`, `.dearmeta*/` 및 각종 임시/리포트 파일이 제외됩니다. 필요 시 CLI를 이용해 언제든지 다시 생성할 수 있습니다.

## 개발 가이드
- 코드 위치: `src/dearmeta/`
- R 분석 스크립트: `scripts/analysis.R`
- 추가 유틸리티 및 설정 예시는 `scripts/`와 `runtime/` 하위 파일을 참고하세요.

테스트는 추후 `tests/` 디렉터리에 추가할 예정입니다. 새로운 기능을 구현할 때는 `ruff`, `pytest` 등 개발 의존성을 활용해 품질을 검증해 주세요.

## 라이선스

이 프로젝트는 Apache License 2.0 하에 배포됩니다. 자세한 사항은 `LICENSE` 파일을 참고하세요.
