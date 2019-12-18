using Logging
using Dates
using Formatting

function fmt(level, _module, group, id, file, line)
    return (:blue, format("{:<23}:", Dates.now()), "")
end

const logger = ConsoleLogger(meta_formatter=fmt)
