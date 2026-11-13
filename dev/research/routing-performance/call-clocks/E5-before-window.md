# E5 — часы Core при восьми workers

Все измеренные времена, поправки и выводы ниже — **exploratory**: по одному
диагностическому процессу P1/P2, без A/B/A. Координатор запускал их вне sandbox,
последовательно в общем heavy slot; light work других сессий продолжалась,
load average был примерно 1,4–3,5. Это не квалифицированное timing window.

P2 отделяет около 1,168 с CPU потоков, запускаемых Core, от CPU Python caller.
В P1 одновременно работающие caller threads и ожидания делают сумму Core spans
намного больше elapsed wall. Инструментация существенно меняет саму нагрузку;
эти часы не доказывают serial floor, free-Core endpoint или конечный выигрыш.
Sensitivity slopes пока отсутствуют: согласованное тихое окно ещё не проведено.

Источники: S = /Users/networmix/ws/work/2026-11-core-routing-research.noindex/scratch/clocks/.
Новый reference extension: SHA-256
208487c7aed158fb418791c22e23d09332533595053c0e4810923b31244f7c84.
NetGraph 8a641a1, Core 9cd7642, CPython 3.14.10 с GIL; входы — неизменённые
orient-add/workloads/dcbb96.yml и nsfnet.yml, workers=8. Старые результаты
ORIENT-ADD не подставляются вместо нового reference.

**Analysis wall и CPU, секунды, exploratory.** Это фаза ngraph.run: MSD и
все placement steps P1; три all-pairs анализа P2. Core span окружает исходный
public call, constructor body, getter или scoped release, включая bindings,
GIL reacquisition и внутреннее ожидание. В P1 учтены 140 lifetime Python threads,
а не 140 одновременно работающих threads; P2 имеет один Python caller и native workers.

| Величина | P1 / 8 | P2 / 8 |
| --- | ---: | ---: |
| Analysis elapsed wall | 3,948728 | 4,785017 |
| Global union Core spans на wall-оси | 2,706664 | 0,566050 |
| Сумма Core span wall по caller threads | 18,970862 | 0,566050 |
| Сумма caller CPU внутри spans | 10,837546 | 0,333507 |
| Сумма caller CPU вне spans, raw | 3,424068 | 4,307027 |
| Вне spans после single-thread null correction | 3,089269 | 0,597954 |
| Сумма outside wall, raw | 6,381284 | 4,218967 |
| Outside wall − outside caller CPU, raw | 2,957216 | −0,088060 |
| То же после null correction | 2,964937 | −0,001563 |
| Весь caller lifetime wall − CPU, сумма | 11,090532 | 0,144482 |
| Process CPU − registered Python CPU, endpoint bracket | 0,006753–0,006760 | 1,168092–1,168098 |
| Python self CPU − Mach CPU | 0,000735 | −0,000007 |
| Public call / release events в analysis | 54 506 | 614 480 |

Последние brackets учитывают порядок endpoint reads, а не всю погрешность и не
95% CI. В P1 малый residual включает неучтённые startup/teardown хвосты Python
threads: приписывать эти 6,8 ms потокам Core нельзя. В P2 Python caller один;
import/load/export residual около нуля, joined-native-worker control сохраняет
его CPU после join, а sampled captures видят native evaluator workers.
Поэтому основной residual P2 относится к запущенным Core threads вместе с их
pthread/runtime расходами. Это не чистый solver CPU и не wall time. Per-call
process deltas записаны, но не суммируются через пересекающиеся вызовы.

В P2 116 evaluate_pairs имеют суммарный API wall 0,197183 с и caller CPU
0,054070 с. Native workers работают внутри этих calls, пока caller может ждать
join; process CPU их сохраняет, даже когда threads уже завершены. Оставшиеся
run events в основном getters конфигурации/результатов. Нельзя делить общий
native-worker CPU на caller span wall и называть результат speedup.

**Погрешность и воздействие сборщика.** Single-thread null wrapper стоит примерно
6,31 мкс CPU/call в P1 и 6,19 мкс в P2; inner null span — около 0,165 и 0,161 мкс CPU.
Полная поправка outside вычитает full-wrapper cost минус inner-span cost.
При восьми registered callers плюс main null cost выше: registry snapshots
читают несколько CPU clocks и синхронизируются RLock. Для P1 применение разных
калибровок меняет outside CPU proxy от 1,825155 до 3,089269 с. Это чувствительность
к поправке, не доверительный интервал и не доказанные границы serial работы.
Восьмипоточная поправка к однопоточному Python пути P2 неприменима; отрицательные
значения такой диагностической чувствительности в JSON не являются временем.

Отрицательные raw outside wall−CPU у P2 объясняются clock-order cost, накопленным
на сотнях тысяч вызовов: span wall включает больше чтений часов, чем span CPU.
После поправки остаются небольшие signed остатки; они не обрезаны до нуля.
Outside CPU содержит Python/runtime, retention и recording, а binding work
с удержанием GIL остаётся внутри spans. RLock ожидания сборщика также могут
попасть в time-not-running. Эта величина включает GIL, futures/queues, native
joins и OS scheduling; отдельное GIL waiting не измерено.

Диагностические фазы import/load/prepare/run/export/cleanup P1:
0,042 / 6,900 / 1,023 / 3,949 / 0,665 / 0,022 с; сумма 12,600 с.
P2: 0,027 / 0,035 / 0,004 / 4,785 / 9,824 / 0,031 с; сумма 14,707 с.
Это не полный CLI: начальные Core/NumPy imports, установка/калибровка,
evidence processing и проверки экспорта вне этих clocks. Для масштаба:
ранние plain controls дали P1 run 2,772 с и P2 run/export 0,486/0,997 с.
Они получены в других процессах/условиях и не являются парным speedup comparison.
Рост времени P2 до нескольких секунд в основном иллюстрирует стоимость
1,93 млн instrumented boundaries. Peak RSS diagnostics: P1 645,6 MiB, P2 429,8 MiB.

**Что это говорит об endpoint model.** Ни global span union, ни сумма native
CPU не являются устранимым wall-бюджетом: они содержат ожидания и перекрываются
с outside работой других callers. Даже исправленный outside CPU лишь proxy
в инструментированном расписании. Условия STATE §3 — serial, fixed outside
work на critical path — не доказаны. Прежние условные потолки P1/P3 остаются
условными; новые часы не превращают их в endpoint measurement. Конечный
multi-worker gain не вычисляется.

**Sensitivity — ожидается окно координатора.** По RULINGS.md δ=0/0,1/0,2 ms
добавляется один раз перед public call только в P1 case-execution threads.
Destructors/scoped releases не задерживаются; load/prepare/MSD/export и P2/P3
не получают δ. Нулевой control проходит тот же GIL-free helper. План: семь
A/A/A triples, затем семь A/B/A на каждую ненулевую δ, fresh processes,
median B/mean(A1,A2), 95% whole-triple bootstrap, seed 20261003. Дополнительно
отдельные secant slopes для analysis, parallel schedule и whole-process wall.
A/A response floor и пригодность окна ещё неизвестны. Pilot использован только
для настройки и бюджета; slope, CI или free-Core inference из него здесь нет.

**Проверка и воспроизведение.** Все пять экспортов координатора, включая clocks
и captures, повторно побайтово сопоставлены с соответствующими plain exports:
равны без исключений metadata, поскольку build и worker count совпадают.
Boundary controls ранее отвергли poisoned values/types/keys/order; реальные
poisoned raw-event и sample-hash controls тоже отказали. Closure проверена
независимыми union/complement расчётами для каждого thread/phase.

| Данные / код | Путь относительно S или task |
| --- | --- |
| Raw clocks и exports | S/measured-p1-clocks-01/, S/measured-p2-clocks-01/; events.i64, python-vectors.i64, summary.json, results.json |
| Анализы с worker/build metadata | S/measured-p1-analysis-02.json, S/measured-p2-analysis-02.json; численные поля идентичны -01 |
| Calibrations | S/measured-null-8-01.json; calibration в raw summaries; S/mach-probe-02.json |
| Пересчитанные таблицы и byte equality | S/report-data-01.json; S/metadata-analysis-rerun-01.json |
| Методы / точные команды | run_workload.py, collector.py, mach_clock.py, analyze_clocks.py, sensitivity.py; HANDOFF.md |
| Условия / полные launch logs | /Users/networmix/ws/work/2026-11-core-routing-research.noindex/scratch/coord-checks/clocks/run-measured-01.log |

Первые новые one-worker P1 clocks ещё не запускались. По RULINGS.md, пункт 7,
N для E1 — global_span_union_all_phases_s: сумма unions всех инструментированных
фаз, включая MSD, getters, scoped releases и неизменённые Core calls.
Analysis-only union в phases.run — отдельный subset; он не заменяет полный N.
Непересечение фаз проверяется assertion. Точная команда и JSON keys — в EXIT.md.
