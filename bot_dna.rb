# frozen_string_literal: true

require 'telegram/bot' # library for telegram API
require 'bio' # library for traslation mRNA and other functions

def cleaner_string(seq) #function of cleaninig the string
  seq = seq.downcase.split('')
  seq.delete_if { |c| c != 'a' and c != 't' and c != 'g' and c != 'c' }
  seq.join('')
end
def get_random_sequence(key)#random function
  random = Random.new
  charset_DNA = ['a', 't', 'g', 'c']
  charset_RNA = ['a', 'u', 'g', 'c']
  charset_PROTEIN = ['Val', 'Asn', 'Thr', 'Ala', 'Leu', 'Cys', 'Asp', 'Glu',  'Phe', 'Gly', 'His', 'Ile', 'Lys', 'Met', 'Pro', 'Gln', 'Arg', 'Ser', 'Trp', 
  'Tyr']
  case key
  when 'ДНК'
  Array.new(random.rand(30...50)) {charset_DNA.sample}.join
  when 'РНК'
  Array.new(random.rand(30...50)) {charset_RNA.sample}.join
  when 'Белок'
  Array.new(random.rand(30...50)) {charset_PROTEIN.sample}.join("-")
  end
end
def dna_to_rna(dna) # my own function that parse the string and make transcription 3' 5'DNA ---> mRNA
  rna = []
  dna = dna.split('')
  i = 0
  while i < dna.length # loop for changing letters
    case dna[i]
    when 'a'
      rna << 'u'
    when 't'
      rna << 'a'
    when 'g'
      rna << 'c'
    when 'c'
      rna << 'g'
    end
    i += 1
  end
  rna.join('') # return new string mRNA
end

def complementation(seq)#DNA - DNA
  seq = Bio::Sequence::NA.new(cleaner_string(seq.downcase))
  seq.complement.to_s.reverse
end

def get_open_read_rna(rna) # i this function i parse mRNA to get open read fragments (START-codone ---- STOP-codone)
  i = 0
  j = 0

  begin_of_translation = 0

  while i < rna.length
    if (rna[i] == 'a') && (rna[i + 1] == 'u') && (rna[i + 2] == 'g') # first, finding AUG codone
      begin_of_translation = i
      break
    end
    i += 1
  end

  rna_codones = rna[begin_of_translation..rna.length].scan(/.../) # then i can use method scan() to turn my sequence like "aatggtat" into ["aat","ggt"]
  end_of_translation = rna_codones.length # this is necessary, if in sequence there is no stop-codones. In that case transltion of protein contunues till the end

  while j < rna_codones.length
    if (rna_codones[j] == 'uaa') || (rna_codones[j] == 'uag') || (rna_codones[j] == 'uga') # finding STOP-codones
      end_of_translation = j
      break
    end
    j += 1
  end

  rna_translate_codones = rna_codones[0..end_of_translation - 1].join('') # making the sting[start..stop]

  rna_translate_codones # return it
end
def translate_from_rna(rna)
rna_open_read = get_open_read_rna(rna)
rna_new= Bio::Sequence::NA.new(rna_open_read)
rna_new.translate.codes.join("-")
end
token = '1352839940:AAErfrfImdrECOdwcyJhLzFXpTSdggVfHGE'

Telegram::Bot::Client.run(token) do |bot|
  bot.listen do |message|
    case message.text
    when '/start'
      bot.api.sendMessage(chat_id: message.chat.id, text: "Здравствуй, #{message.from.first_name}.\nЯ могу анализировать последовательности ДНК.\nЧтобы получить список команд, введите /help")
    when '/help'
      bot.api.sendMessage(chat_id: message.chat.id, text: "Вот мои команды:
/full_analysis_5_3 - делает полный анализ смысловой(5' 3') последовательности ДНК (5'ДНК 3' - 3'ДНК 5'- иРНК - Белок) \n /full_analysis_3_5 - делает полный анализ транскрибируемой(3' 5') последовательности ДНК (3'ДНК 5'- иРНК - Белок)\n /get_random - создает рандомную последовательность выбранного типа\n /info - Дополнительная информация, которая может быть полезной\n /complement - выводит цепь, которая комплиментарна введённой \n /translate - проводит трансляцию введенной последовательности иРНК")

    when '/full_analysis_5_3'
      bot.api.sendMessage(chat_id: message.chat.id, text: 'Вводите последовательность ДНК')
      bot.listen do |dna|
        _dna = cleaner_string(dna.text)
        if _dna.length < 5
          bot.api.sendMessage(chat_id: message.chat.id, text: 'Ваша последовательность слишком короткая или в ней были допущены ошибки')
          break
        else
          DNA = Bio::Sequence::NA.new(_dna)
          dna_3_5 = DNA.complement.to_s.reverse.upcase.split('').join(' - ')
          dna_5_3 = DNA.to_s.upcase.split('').join(' - ')
          rna = dna_to_rna(DNA.complement.to_s.reverse.downcase).to_s
          protein = Bio::Sequence::NA.new(get_open_read_rna(rna).to_s)
          protein = protein.translate.codes.join('-')

          File.open('result.txt', 'w') { |f| f.write "Смысловая цепь ДНК: \n5' #{dna_5_3} 3'\nТранскрибируемая цепь ДНК: \n 3' #{dna_3_5} 5' \n ---> \n иРНК: \n5' #{rna.upcase.split('').join(' - ')} 3' \n Кодирующие триплеты иРНК:\n---> \n #{get_open_read_rna(rna).scan(/.../).join('-').upcase}\n ---> \n Последовательность анимокислот в белке: \n #{protein}" }
          bot.api.sendMessage(chat_id: message.chat.id, text: 'Вот результаты исследования')
          bot.api.sendDocument(chat_id: message.chat.id, document: Faraday::UploadIO.new('result.txt', 'text/plain'))
          File.open('result.txt', 'w') { |file| file.truncate(0) }

          break
        end
      end

    when '/full_analysis_3_5'
      bot.api.sendMessage(chat_id: message.chat.id, text: 'Введите транскрибируемую последовательность ДНК')

      bot.listen do |dna|
        if cleaner_string(dna.text).length < 5
          break
        else
          DNA = Bio::Sequence::NA.new(cleaner_string(dna.text))
          dna_3_5 = DNA.to_s.upcase.split('').join(' - ')
          rna = dna_to_rna(DNA.to_s.downcase).to_s
          protein = Bio::Sequence::NA.new(get_open_read_rna(rna).to_s)
          protein = protein.translate.codes.join('-')

          File.open('result.txt', 'w') { |f| f.write "\nТранскрибируемая цепь ДНК: \n 3' #{dna_3_5} 5' \n ---> \n иРНК: \n5' #{rna.upcase.split('').join(' - ')} 3' \n Кодирующие триплеты иРНК:\n---> \n #{get_open_read_rna(rna).scan(/.../).join('-').upcase} \n---> \n Последовательность анимокислот в белке: \n #{protein}" }
          bot.api.sendMessage(chat_id: message.chat.id, text: 'Вот результаты исследования')
          bot.api.sendDocument(chat_id: message.chat.id, document: Faraday::UploadIO.new('result.txt', 'text/plain'))
          File.open('result.txt', 'w') { |file| file.truncate(0) }

          break

        end
      end
    when '/info'
     bot.api.sendDocument(chat_id: message.chat.id, document: Faraday::UploadIO.new('info.txt', 'text/plain'))
  when '/get_random'
      question = 'Введите тип последовательности(Белок, ДНК, РНК)'
    # See more: https://core.telegram.org/bots/api#replykeyboardmarkup
    answers = Telegram::Bot::Types::ReplyKeyboardMarkup.new(keyboard: ['Белок', 'ДНК', 'РНК'], one_time_keyboard: true)
    bot.api.sendMessage(chat_id: message.chat.id, text: question, reply_markup: answers)
    bot.listen do |answers|
      case answers.text
  when 'Белок'
    bot.api.sendMessage(chat_id: message.chat.id, text: get_random_sequence('Белок'))
  when 'ДНК'
    bot.api.sendMessage(chat_id: message.chat.id, text: get_random_sequence('ДНК'))
  when 'РНК'
    bot.api.sendMessage(chat_id: message.chat.id, text: get_random_sequence('РНК'))
  else 
    break
  end
      break
    end

when '/complement'
  bot.api.sendMessage(chat_id: message.chat.id, text: "Введите последовательность")
  bot.listen do |seq|
  bot.api.sendMessage(chat_id: message.chat.id, text: complementation(seq.text))
break
end
when '/translate'
  bot.api.sendMessage(chat_id: message.chat.id, text: "Введите последовательность")
  bot.listen do |seq|
  bot.api.sendMessage(chat_id: message.chat.id, text: translate_from_rna(seq.text))
break
end
end
  end
end

